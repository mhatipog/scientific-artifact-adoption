# Round 8: first direct microscope-image audit

One Q1 water–air TIFF file from the CC-BY DaRUS dataset (Yan et al.
2025, doi:10.18419/DARUS-4683) is selected by public metadata, streamed
into ephemeral GitHub Actions storage, checked against the archive MD5,
and sampled at twelve frames. A derived contact sheet and small JSON
summary are stored as short-retention workflow artifacts. The original
video bytes are **not committed to GitHub or republished**.

The archive names sequences using **original video frame numbers**;
e.g. "60_Q1_water_air_3330-3710_TD.tif". The paper states 30 fps
and archive metadata describes ~2.4859 µm/pixel. The corner radius
5.6 µm corresponds to only 2.25 pixels, and the earlier planar
2.93 µm reversal-distance estimate to only 1.18 pixels, making
subpixel measurement and an uncertainty budget essential.

**Do not infer synchronized camera–pressure event times** simply from
the frame filename: the pressure logger has an independent absolute
start clock, and the TIFF may lack acquisition timestamps. Video
clock alignment, real advancing angles, actual interface segmentation,
3D geometry and apparatus compliance remain outstanding. A contact sheet
alone is not a result validating our theory.

Workflow is manual and PR-path-scoped; no schedule. Avoid expensive
large-data downloads in every pull request unrelated to CAH imagery.
