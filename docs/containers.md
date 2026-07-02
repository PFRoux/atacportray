# Custom containers

Most processes in atacportray use public biocontainers pulled automatically by
Nextflow. Three tools have no suitable public image and ship a `Dockerfile`
next to their module:

| Module | Image tag | Why custom |
|--------|-----------|------------|
| ROSE | `ghcr.io/stjude/abralab/rose:v1.3.2` | Uses the public St. Jude ABRALab image (no build needed) |
| TelomereHunter | `ghcr.io/pfroux/atacportray-telomerehunter:1.1.0` | Python 2.7 only, not on bioconda |
| NucleoATAC | `ghcr.io/pfroux/atacportray-nucleoatac:1.0.0` | Upstream is Python 2.7; built from the Python 3.11 community fork (PR #99) |

ROSE needs no action — it pulls a pre-built public image. **TelomereHunter and
NucleoATAC must be built once and pushed to GHCR** before you run their branches
with `-profile docker` (or generate Singularity images from those tags).

## One-time build & push

You need a GitHub Personal Access Token with `write:packages` scope. Log in to
GHCR:

```bash
echo "$GHCR_PAT" | docker login ghcr.io -u PFRoux --password-stdin
```

Build and push both images (run from the repo root):

```bash
# TelomereHunter (python2.7 + R plotting)
docker build --platform=linux/amd64 \
    -t ghcr.io/pfroux/atacportray-telomerehunter:1.1.0 \
    modules/local/telomerehunter
docker push ghcr.io/pfroux/atacportray-telomerehunter:1.1.0

# NucleoATAC (python3.11 fork, GreenleafLab PR #99)
docker build --platform=linux/amd64 \
    -t ghcr.io/pfroux/atacportray-nucleoatac:1.0.0 \
    modules/local/nucleoatac
docker push ghcr.io/pfroux/atacportray-nucleoatac:1.0.0
```

After the first push, make each package **public** in your GitHub account
(Packages → the package → Package settings → Change visibility → Public) so
Nextflow can pull it without credentials on any machine.

> On Apple Silicon the `--platform=linux/amd64` flag is required (both tools
> depend on amd64-only conda packages); Docker runs them under emulation.

## Singularity / Apptainer

The module container directives already provide `oras://ghcr.io/...` addresses
for Singularity. Once the images are public on GHCR, `-profile singularity`
pulls and converts them automatically — no separate build step.
