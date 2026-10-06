#!/usr/bin/env bash
# Rerun every command of the documentation's worked examples (docs/examples/).
#
# Usage: bash examples/run_worked_examples.sh [OUTPUT_DIR]   (default results/worked-examples)
#
# Needs CREVICE installed (the `crevice` command), network access to download
# 1GRM, 1OED, 3UKM (biological assembly 1) and 4PYP from RCSB, about 15
# minutes and about 50 MB. Each example runs in its own subdirectory and
# appends the commands and their output to OUTPUT_DIR/commands.log; the
# structures are cached once in OUTPUT_DIR/.crevice/pdb. Viewer images are not
# rendered (open the .pml/.vmd/.cxc scenes in your viewer).
set -euo pipefail

out="$(mkdir -p "${1:-results/worked-examples}" && cd "${1:-results/worked-examples}" && pwd)"
cache="$out/.crevice/pdb"
log="$out/commands.log"
: > "$log"

run() {  # run DIR crevice-arguments...
    local dir="$1"; shift
    mkdir -p "$out/$dir"
    echo "[$dir] \$ crevice $*" | tee -a "$log"
    (cd "$out/$dir" && crevice "$@") 2>&1 | tee -a "$log"
}

crevice fetch 1GRM 1OED 4PYP --cache-dir "$cache" | tee -a "$log"
crevice fetch 3UKM --assembly 1 --cache-dir "$cache" | tee -a "$log"
crevice fetch 1GRM --format pdb --cache-dir "$cache" | tee -a "$log"

# 1GRM: pore radius profile and channel cast
run 1grm profile 1GRM --cache-dir "$cache" --axis auto -o 1GRM_profile.json \
    --png 1GRM_profile.png --annotated-png 1GRM_profile_residues.png --dpi 150
run 1grm publish 1GRM --cache-dir "$cache" --out-dir 1GRM_bundle --prefix 1GRM --skip-hydration --dpi 150

# 3UKM: a capped channel with lateral exits
run 3ukm profile 3UKM --assembly 1 --cache-dir "$cache" --exclude-hetero --lateral-exits \
    -o 3UKM_profile.json --png 3UKM_profile.png --dpi 150
run 3ukm publish 3UKM --assembly 1 --cache-dir "$cache" --exclude-hetero --lateral-exits \
    --enclosure-radius 0.55 --skip-hydration --skip-cavities --skip-tunnels --skip-network \
    --out-dir 3UKM_bundle --prefix 3UKM --dpi 150

# 1OED: a wide pore with gaps in its wall
run 1oed profile 1OED --cache-dir "$cache" -o 1OED_profile.json --png 1OED_profile.png --dpi 150
run 1oed profile 1OED --cache-dir "$cache" --enclosure-radius 2.0 -o 1OED_profile_2.0.json
run 1oed residues 1OED --cache-dir "$cache" --enclosure-radius 2.0 -o 1OED_residues.csv \
    --png 1OED_residues.png --dpi 150

# 4PYP: interior cavity cast and its boundary residues (the profile is expected to fail)
run 4pyp profile 4PYP --cache-dir "$cache" -o 4PYP_profile.json || echo "[4pyp] profile unresolved, as expected" | tee -a "$log"
run 4pyp cast 4PYP --cache-dir "$cache" --out-dir 4PYP_cast --skip-hydration
run 4pyp residue-evidence 4PYP_cast/4PYP_viewer.pdb --volume-dx 4PYP_cast/4PYP_volume.dx \
    --scene 4PYP_cast/4PYP_volume.pml --out-dir 4PYP_evidence --prefix 4PYP --skip-hydration --dpi 150

# Atomic radii and figure text
mkdir -p "$out/radii"
cat > "$out/radii/my_radii.json" <<'JSON'
{
  "name": "bondi-larger-oxygen",
  "base": "bondi",
  "element_radii": {"O": 1.60},
  "citation": "worked example: O 1.52 -> 1.60 A to show the effect of one radius"
}
JSON
run radii profile 1GRM --cache-dir "$cache" --axis auto --enclosure-radius 0.75 -o standard.json --png standard.png --dpi 150
run radii profile 1GRM --cache-dir "$cache" --axis auto --enclosure-radius 0.75 --radii hole -o hole.json
run radii profile 1GRM --cache-dir "$cache" --axis auto --enclosure-radius 0.75 --radii my_radii.json -o custom.json
run radii profile 1GRM --cache-dir "$cache" --axis auto --enclosure-radius 0.75 --annotate -o annotated.json \
    --png annotated.png --dpi 150

# A small public trajectory: the 1GRM NMR ensemble
run trajectory trajectory --topology "$cache/1GRM.pdb" --selection all --inspect "$cache/1GRM.pdb"
run trajectory trajectory "$cache/1GRM.pdb" --topology "$cache/1GRM.pdb" --selection all --skip-hydration \
    -o 1GRM_nmr.json --png 1GRM_nmr_profiles.png --distribution-csv 1GRM_nmr_distribution.csv \
    --network-json 1GRM_nmr_network.json --dpi 150

echo "Done. Outputs in $out; commands and output in $log"
