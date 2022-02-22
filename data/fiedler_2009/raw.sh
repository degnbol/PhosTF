#!/usr/bin/env zsh
mkdir -p raw
cd raw
# see https://www.cell.com/supplemental/S0092-8674(09)00002-6#supplementaryMaterial
# Table S1. All Genes in E-MAP
wget 'https://www.cell.com/cms/10.1016/j.cell.2008.12.039/attachment/479ba893-11ad-4229-af0f-f3ac149d5387/mmc2.xls'
# Table S2. E-MAP Scores
wget 'https://www.cell.com/cms/10.1016/j.cell.2008.12.039/attachment/cc2b6a4d-abb1-41ce-86cf-fd812eab882d/mmc3.txt'
# Table S3. Known Phosphorylation/Dephosphorylation Events, Assembled from the Literature
wget 'https://www.cell.com/cms/10.1016/j.cell.2008.12.039/attachment/f795bf08-10c0-412d-8a69-55b63e13efd4/mmc4.xls'
# Table S4. Scores for Known Phosphorylation Events
wget 'https://www.cell.com/cms/10.1016/j.cell.2008.12.039/attachment/5a8815ca-46ac-4273-918f-533d6442ed18/mmc5.xls'
# Table S5. Ratios for Figure 2B
wget 'https://www.cell.com/cms/10.1016/j.cell.2008.12.039/attachment/9089f0d9-18bb-4b12-8138-e7786a9a7a99/mmc6.xls'
# Table S6. Most Positive Kinase-Phosphatase Pairs
wget 'https://www.cell.com/cms/10.1016/j.cell.2008.12.039/attachment/3d692020-5d8d-4f6f-937c-e70822701e45/mmc7.xls'
# Table S7. TGMs
wget 'https://www.cell.com/cms/10.1016/j.cell.2008.12.039/attachment/80b3198b-c87a-47de-be24-60da2b93829f/mmc8.xls'
# Table S8. Gene Expression Data
wget 'https://www.cell.com/cms/10.1016/j.cell.2008.12.039/attachment/e204439e-2d6a-41fc-a10a-58b49d82f80a/mmc9.xls'
cd -
