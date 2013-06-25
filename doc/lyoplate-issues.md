# Lyoplate 3.0 Issues

### Problems with Samples

* Miami did not include FSC-H in their FCS files. To standardize the gating
  schemes across all centers, we did not apply a singlet gate. This had little
  effect on the downstream gates and the corresponding proportions.

* FCS files from CIMR had a varying number of markers for the T-cell panel. To
  handle this case, we had to semi-automatically standardize the FCS files.

* We removed BSMS from all panels because they failed to include key markers for
  four panels:
  1. B-cell: IgD
  2. T-reg: CD25
  3. DC/Mono/NK: CD56
  4. T-helper: CXCR3

* The Lyoplate protocol mentions that the FSC threshold should be kept low
  enough so that some of the debris is visible. The interpretation of this
  statement varied widely. For example, approximately 55% of the cells in sample
  12828 from Yale were debris. Contrarily, hardly any debris were present in
  Miami. Despite the varying amounts of debris, we were able to gate the
  lymphocyte populations without difficulty.

### Compensation and Transformation

* Perhaps the most challenging aspect of the automated analysis is to determine
  proper transformation.  We applied several standard transformation methods to
  the samples and observed poor transformations in many cases. We found that
  pregating the samples often improved the transformations. However, even with
  the removal of debris and boundary events, the transformations proved to be
  poor for several markers. In several cases (e.g, CD27/IgD), our automated
  gating pipelines were robust against poor transformation. Contrarily, there
  were instances where the poor transformation yielded data that were quite
  challenging to gate (e.g, CD4/CXCR3 vs CD4/CCR6 in the T-helper panel). We
  would have been ideally provided with the same transformations as those used
  in the manual analysis. This would remove a source of variability between the
  manual and automated gating. Furthermore, in several instances the gating
  would have been less challenging and more straightforward.

* The manner in which the compensation controls were provided was nonstandard.
  - For some centers, such as Stanford, generic color controls were used for
  compensation. In some of these cases, a single FCS file was used for all
  markers sharing the same channel.
  - For other centers, such as NHLBI, one compensation control was used for
  several markers that share the same channel (e.g., NHLBI), while other centers
  use one compensation control per marker (e.g., CIMR). To handle this, we use
  only the first FCS file if the channel names were non-unique.
  - BSMS did not provide an unstained compensation control. In this case, the
  single-stained controls doubled as unstained controls because the stained and
  unstained peaks were well-defined. We extracted the negative peaks from one of
  the controls and use it as an unstained control. No mention was given
  initially to indicate that this is how an unstained compensation control
  should be extracted. Only after we inquired did BSMS provide the
  aforementioned strategy for extracting an unstained control.

* The spillover matrices for Yale are likely unreliable. For the FITC-A
  compensation control, there were approximately 137 cells actually stained.
  Furthermore, for several channels within this control, there appeared to be a
  small number of stained cells. This likely explains the poor spillover
  matrix. The spillover via Yale's PerCP-Cy5-5-A compensation control is likely
  unreliable too, although not as extreme as FITC-A.

### Excel Configuration Files

* Several cases where the filenames provided in the Excel files did not match
  the actual filenames. In these cases, the Excel files were edited manually to
  ensure the files matched correctly.

* Channel and marker names were often not uniform nor did they match those given
  in the FCS files. We processed FCS files and semi-automatically converted the
  channel and marker names. For example, we used the name Lineage to standardize
  CD3CD19CD20, CD3+19+20, CD3_CD19_CD20, CD3+CD19+CD20+, CD3+CD19+CD20, and
  CD3+19+20.

* Related to the previous point, several marker and channel names were not easily
  parsed. We manually edited these entries to automate the downstream
  preprocessing and gating.
