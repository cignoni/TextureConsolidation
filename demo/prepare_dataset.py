"""
Convert a dataset (textured PLY + cameras.json + images) into the demo layout:

    <out>/mesh.ply          binary PLY, per-face vertex_indices + texcoord
    <out>/cameras.json      frames in view order, file_path = images/<i>.png
    <out>/images/<i>.png    i = 0 .. N-1, in cameras.json order

and repair the mesh so TextureDomainStitching accepts it:
  * drop faces with zero area in texture space;
  * split "pinch" vertices, i.e. vertices whose incident faces form more than one
    edge-connected fan (TSP aborts with "Bad components" on these);
  * remove vertices no face uses.
With --reparameterize, the texture coordinates are replaced by a fresh xatlas atlas (the 3D vertices and
triangles are kept), for meshes whose original atlas TextureSignalProcessing cannot handle.
With --area A, the mesh is uniformly scaled to surface area A (camera translations are scaled to match, so the
images are unchanged). The flow parameters depend on mesh scale; the release defaults were tuned at area ~5.
Images can be downsampled by an integer factor; intrinsics are scaled to match.

Usage: python demo/prepare_dataset.py <in_dir> <out_dir> [--mesh in.ply] [--downsample 2] [--reparameterize] [--area 5]
Needs numpy and opencv-python (and xatlas for --reparameterize).
"""
import argparse, json, shutil
from pathlib import Path
import numpy as np
import cv2


def read_ply(path):
    data = path.read_bytes()
    end = data.index(b"end_header\n") + len(b"end_header\n")
    header = data[:end].decode()
    if "format binary_little_endian" not in header:
        raise RuntimeError("expected a binary_little_endian PLY")
    lines = header.splitlines()
    nv = int(next(l for l in lines if l.startswith("element vertex")).split()[-1])
    nf = int(next(l for l in lines if l.startswith("element face")).split()[-1])
    vtype = "f4" if "property float x" in header else "f8"
    vprops = [l for l in lines if l.startswith("property") and lines.index(l) < lines.index(next(l for l in lines if l.startswith("element face")))]
    if len(vprops) != 3:
        raise RuntimeError(f"expected vertex properties x,y,z only, got {vprops}")
    vbytes = nv * 3 * np.dtype(vtype).itemsize
    V = np.frombuffer(data[end:end + vbytes], dtype="<" + vtype).reshape(nv, 3).astype(np.float64)
    fdt = np.dtype([("n", "u1"), ("v", "<i4", 3), ("m", "u1"), ("uv", "<f4", 6)])
    F = np.frombuffer(data[end + vbytes:end + vbytes + nf * fdt.itemsize], dtype=fdt)
    if len(F) != nf or (F["n"] != 3).any() or (F["m"] != 6).any():
        raise RuntimeError("unexpected face layout (expected list uchar int vertex_indices, list uchar float texcoord with 6 entries)")
    comments = [l for l in lines if l.startswith("comment")]
    return V, F["v"].copy(), F["uv"].reshape(-1, 3, 2).copy(), comments


def write_ply(path, V, T, UV, comments):
    fdt = np.dtype([("n", "u1"), ("v", "<i4", 3), ("m", "u1"), ("uv", "<f4", 6)])
    F = np.zeros(len(T), dtype=fdt)
    F["n"], F["v"], F["m"], F["uv"] = 3, T, 6, UV.reshape(-1, 6)
    header = ["ply", "format binary_little_endian 1.0", *comments,
              f"element vertex {len(V)}", "property float x", "property float y", "property float z",
              f"element face {len(T)}", "property list uchar int vertex_indices", "property list uchar float texcoord",
              "end_header"]
    with open(path, "wb") as f:
        f.write(("\n".join(header) + "\n").encode())
        f.write(V.astype("<f4").tobytes())
        f.write(F.tobytes())


def split_pinch_vertices(V, T, UV):
    """Give each edge-connected fan of faces around a texture-space vertex its own copy of the vertex.

    TextureSignalProcessing identifies a texture-space vertex by (3D vertex, texture coordinate), and two
    faces are adjacent in texture space only if they share two such vertices. A (vertex, uv) point whose
    faces form more than one adjacent fan is a pinch that it rejects ("Bad components").
    """
    key = lambda t, k: (int(T[t][k]), float(UV[t][k][0]), float(UV[t][k][1]))
    incident = {}
    for t in range(len(T)):
        for k in range(3):
            incident.setdefault(key(t, k), []).append((t, k))
    V = list(V)
    T = T.copy()
    split = 0
    for faces in incident.values():
        if len(faces) < 2:
            continue
        # union faces that share a texture-space edge through this point (i.e. share another texture-space vertex)
        parent = list(range(len(faces)))
        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]; i = parent[i]
            return i
        owner = {}
        for i, (t, k) in enumerate(faces):
            for other in (key(t, (k + 1) % 3), key(t, (k + 2) % 3)):
                if other in owner:
                    parent[find(i)] = find(owner[other])
                else:
                    owner[other] = i
        fans = {}
        for i in range(len(faces)):
            fans.setdefault(find(i), []).append(faces[i])
        for fan in list(fans.values())[1:]:
            v = T[faces[0][0]][faces[0][1]]  # the original vertex (the first fan keeps it)
            V.append(V[v]); nv = len(V) - 1; split += 1
            for t, k in fan:
                T[t][k] = nv
    return np.array(V), T, split


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("in_dir"); ap.add_argument("out_dir")
    ap.add_argument("--mesh", default="in.ply")
    ap.add_argument("--downsample", type=int, default=1)
    ap.add_argument("--reparameterize", action="store_true", help="replace the texture atlas with an xatlas one")
    ap.add_argument("--area", type=float, default=0, help="scale the mesh (and cameras) to this surface area")
    a = ap.parse_args()
    src, dst = Path(a.in_dir), Path(a.out_dir)
    (dst / "images").mkdir(parents=True, exist_ok=True)

    V, T, UV, comments = read_ply(src / a.mesh)
    if a.reparameterize:
        import xatlas
        _, indices, uvs = xatlas.parametrize(V.astype(np.float32), T.astype(np.uint32))
        UV = uvs[indices].astype(np.float32)   # xatlas keeps the face order, so per-corner uvs map onto T
        print(f"reparameterized: {len(np.unique(indices))} atlas vertices")
    e1, e2 = UV[:, 1] - UV[:, 0], UV[:, 2] - UV[:, 0]
    keep = np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]) > 0
    T, UV = T[keep], UV[keep]
    V, T, split = split_pinch_vertices(V, T, UV)
    used = np.zeros(len(V), bool); used[T.ravel()] = True
    remap = -np.ones(len(V), int); remap[used] = np.arange(used.sum())
    V, T = V[used], remap[T]
    scale = 1.
    if a.area > 0:
        area = 0.5 * np.linalg.norm(np.cross(V[T[:, 1]] - V[T[:, 0]], V[T[:, 2]] - V[T[:, 0]]), axis=1).sum()
        scale = np.sqrt(a.area / area)
        V = V * scale
        print(f"scaled mesh by {scale:.4f}: surface area {area:.4f} -> {a.area}")
    comments = [c for c in comments if not c.startswith("comment TextureFile")] + ["comment TextureFile texture.png"]
    write_ply(dst / "mesh.ply", V, T, UV, comments)
    print(f"mesh: dropped {int((~keep).sum())} zero-UV faces, split {split} pinch vertices, "
          f"removed {int((~used).sum())} unused vertices -> {len(V)} vertices / {len(T)} faces")

    frames = json.loads((src / "cameras.json").read_text())["frames"]
    s = a.downsample
    for i, fr in enumerate(frames):
        im = cv2.imread(str(src / fr["file_path"]), cv2.IMREAD_UNCHANGED)
        if s > 1:
            im = cv2.resize(im, (im.shape[1] // s, im.shape[0] // s), interpolation=cv2.INTER_AREA)
            for key in ("fx", "fy", "cx", "cy"):
                fr[key] = fr[key] / s
            fr["w"], fr["h"] = im.shape[1], im.shape[0]
        cv2.imwrite(str(dst / "images" / f"{i}.png"), im)
        fr["file_path"] = f"images/{i}.png"
        for r in range(3): fr["w2c"][r][3] *= scale
    # one key per line and brackets on their own tokens, as the cameras.json reader expects
    keys = ["w", "h", "fx", "fy", "cx", "cy", "w2c", "file_path"]
    frames = [{k: fr[k] for k in keys} for fr in frames]
    (dst / "cameras.json").write_text(json.dumps({"frames": frames}, indent=4) + "\n")
    print(f"cameras: {len(frames)} views, images {frames[0]['w']}x{frames[0]['h']}")


if __name__ == "__main__":
    main()
