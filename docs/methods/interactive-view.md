# Interactive view surfaces

This note describes what [`crevice view`](../tools/view.md) and
{func}`crevice.view` draw, and how it relates to the native PyMOL, VMD and
ChimeraX scenes. The view only displays results; it measures nothing.

## Surfaces

Every cast surface in an interactive view is the **native display surface**:
the 0.5 isosurface of the same display field CREVICE writes for the native
viewers (`PREFIX_<segment>_display.dx` for ENTRY, LUMEN and EXIT segments,
`PREFIX_display.dx` for a single cast). The mesh is extracted with
marching cubes (`skimage.measure.marching_cubes`) at level 0.5, in the
field's world coordinates. Because the field already includes the run's
display smoothing ([Display surface smoothing](../SURFACE_SMOOTHING.md),
`--smooth`), the browser shows the same surface as PyMOL, VMD and ChimeraX.

For results passed in memory (`crevice.view(frame=..., cast=...)`), the
display field is built from the cast's binary grid with the same smoothing
rules (`smooth=`), on the cast's own grid and without the short display
continuation that `publish` adds beyond each mouth.

Measured volumes always come from the binary `*_volume.dx` grids and are
reported in the tables. A display surface is never used for a measurement.

## Other layers

- **Protein:** a cartoon, grey-blue at 25% transparency.
- **Segments:** ENTRY yellow, LUMEN teal and EXIT rust; a cast that is not a
  channel is violet, as in the native scenes.
- **Lining residues** (hidden at first): sticks with carbons in their
  segment's colour, from the bundle's per-segment lining table.
- **Profile centre line** (hidden at first): tubes through the axial profile
  samples and any side-exit legs.
- **Region outputs** (`region-trajectory`, `cavity-trajectory`, `hydration`):
  the region cast at 35% opacity and the mean water density surface at 60%
  opacity, more opaque than in the native scenes (18% and 42%) so they stay
  visible through the cartoon in a browser.
- **Member waters:** orange spheres for the water oxygens inside the measured
  region in one saved frame (first, middle or last; `--snapshot`). Membership
  is geometric, as defined in
  [Residue hydration and solvent exposure](../HYDRATION_METHODS.md); it is not
  evidence of binding.

The starting orientation matches the native scenes (channel axis vertical).
Nothing is drawn into the 3D picture: layer names appear only as checkbox
labels beside it.

## Rendering

Both renderers use one scene description and one loader script on top of
3Dmol.js (Rego & Koes 2015). The HTML file inlines a vendored 3Dmol.js 2.5.5
build (BSD-3-Clause, licence text embedded), so it opens offline; `--cdn`
loads the library from jsDelivr instead. In Jupyter, `.show()` uses py3Dmol
(the `view` extra) and loads 3Dmol.js from a CDN unless `offline=True`.

The page decodes its compressed meshes with the browser's
`DecompressionStream`, available in Firefox 113, Chrome 80, Safari 16.4 and
later. The software tests check the scene contents and the HTML file; how the
scene looks was checked from screenshots in headless Firefox, and other
browsers were not tried.

## References

- Rego N, Koes D (2015) 3Dmol.js: molecular visualization with WebGL.
  *Bioinformatics* 31:1322–1324.
