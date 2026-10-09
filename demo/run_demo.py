"""
Runs the texture-consolidation pipeline on a dataset:

    1. BackProjection           renderings/photos -> per-view texture and confidence images
    2. TextureFlow              registers the per-view textures with surface optical flow
    3. TextureDomainStitching   fuses the registered textures into one seamless texture

A dataset is a folder containing
    mesh.ply        texture-atlased triangle mesh (binary PLY with per-face "texcoord")
    cameras.json    {"frames": [{"w","h","fx","fy","cx","cy","w2c","file_path"}, ...]}
                    with file_path relative to the dataset folder
and the images it references. Views are numbered 0..N-1 in cameras.json order.

Outputs (in --out):
    backprojection/<i>.png, <i>.confidence.png   per-view back-projected textures
    flow/<i>.png, <i>.confidence.png             registered textures
    texture.png                                  the consolidated texture
    textured_mesh.ply                            the mesh, referencing texture.png (opens in MeshLab)
    blended.png, blended_mesh.ply                (with --baseline) confidence-weighted blend of the
                                                 back-projections, without registration or stitching

Usage:
    python demo/run_demo.py --data <dataset folder> --bin <folder with the executables> [--out demo_output]
"""
import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


def run(cmd, log):
    print("  " + " ".join(str(c) for c in cmd), flush=True)
    with open(log, "w") as f:
        result = subprocess.run([str(c) for c in cmd], stdout=f, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        sys.exit(f"[ERROR] {Path(cmd[0]).name} failed (exit code {result.returncode}); see {log}")


def find_executable(bin_dir, name):
    for candidate in (bin_dir / name, bin_dir / (name + ".exe"), bin_dir / "Release" / (name + ".exe")):
        if candidate.is_file():
            return candidate.resolve()
    sys.exit(f"[ERROR] Could not find {name} in {bin_dir}")


def textured_mesh(mesh, texture, out):
    """Copies a PLY mesh, setting its 'comment TextureFile' to the given texture (relative to the mesh)."""
    data = mesh.read_bytes()
    end = data.index(b"end_header\n")
    lines = [l for l in data[:end].decode().split("\n") if l and not l.startswith("comment TextureFile")]
    lines.insert(2, "comment TextureFile " + os.path.relpath(texture, out.parent).replace("\\", "/"))
    out.write_bytes(("\n".join(lines) + "\n").encode() + data[end:])


def blend(textures, confidences, out):
    """Confidence-weighted average of the back-projected textures (baseline)."""
    import numpy as np
    import cv2
    num = den = 0
    for t, c in zip(textures, confidences):
        w = cv2.imread(str(c), cv2.IMREAD_GRAYSCALE).astype(np.float64)[..., None] / 255.
        num = num + cv2.imread(str(t)).astype(np.float64) * w
        den = den + w
    cv2.imwrite(str(out), np.where(den > 0, num / np.maximum(den, 1e-12), 0).astype(np.uint8))


def main():
    ap = argparse.ArgumentParser(description="Texture consolidation demo")
    ap.add_argument("--data", required=True, type=Path, help="dataset folder (mesh.ply, cameras.json, images)")
    ap.add_argument("--bin", required=True, type=Path, help="folder containing the built executables")
    ap.add_argument("--out", default=Path("demo_output"), type=Path, help="output folder")
    ap.add_argument("--res", type=int, default=1024, help="texture resolution")
    ap.add_argument("--baseline", action="store_true", help="also write a confidence-weighted blend without registration (needs numpy and opencv-python)")
    a = ap.parse_args()

    mesh, cameras = (a.data / "mesh.ply").resolve(), (a.data / "cameras.json").resolve()
    for f in (mesh, cameras):
        if not f.is_file(): sys.exit(f"[ERROR] Missing {f}")
    exe = {name: find_executable(a.bin, name) for name in ("BackProjection", "TextureFlow", "TextureDomainStitching")}

    out = a.out.resolve()
    bp, flow = out / "backprojection", out / "flow"
    for d in (bp, flow):
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True)

    start = time.time()
    print("[1/3] Back-projection")
    run([exe["BackProjection"], "--mesh", mesh, "--cameras", cameras, "--res", a.res, a.res,
         "--outTextures", bp / "%d.png", "--outConfidences", bp / "%d.confidence.png", "--verbose"], out / "backprojection.log")

    print("[2/3] Registration (surface optical flow)")
    run([exe["TextureFlow"], "--mesh", mesh, "--inTextures", bp / "%d.png", "--inConfidences", bp / "%d.confidence.png",
         "--outTextures", flow / "%d.png", "--outConfidences", flow / "%d.confidence.png", "--verbose"], out / "flow.log")

    print("[3/3] Gradient-domain stitching")
    # --multi and --wType 1 are required: they select the confidence-weighted stitching of Eq. 7.
    run([exe["TextureDomainStitching"], "--in", mesh, flow / "%d.png", "--mask", flow / "%d.confidence.png",
         "--multi", "--wType", 1, "--erode", 1, "--erodeEps", 0.1, "--out", out / "texture.png"], out / "stitching.log")
    # TextureDomainStitching can report some errors with exit code 0, so also check that it wrote its output
    if not (out / "texture.png").is_file():
        sys.exit(f"[ERROR] TextureDomainStitching did not write {out / 'texture.png'}; see {out / 'stitching.log'}")
    textured_mesh(mesh, out / "texture.png", out / "textured_mesh.ply")

    if a.baseline:
        n = len(list(bp.glob("*.confidence.png")))
        blend([bp / f"{i}.png" for i in range(n)], [bp / f"{i}.confidence.png" for i in range(n)], out / "blended.png")
        textured_mesh(mesh, out / "blended.png", out / "blended_mesh.ply")

    print(f"Done in {time.time() - start:.0f} s. Open {out / 'textured_mesh.ply'} in MeshLab to view the result.")


if __name__ == "__main__":
    main()
