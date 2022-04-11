# dp_tools

# Very rough guidelines and order of steps
1. Define the components (assay related, but intended to be reusable, e.g. ReadsComponent should fit a number of NGS assays) including what data assets are attached [Components from BulkRNASeq Example](https://github.com/J-81/dp_tools/blob/main/dp_tools/components/components.py)
1. Define the sample (assay specific) including what components should be attached 
1. Define the dataset (assay specific) including what components should be attached (e.g. DifferentialGeneExpressionComponent is a dataset component rather than tied to a specific sample) [BulkRNASeq Entities Example](https://github.com/J-81/dp_tools/blob/main/dp_tools/bulkRNASeq/entity.py)
1. Write a loader that connects file location knowledge and a root directory input to assemble the dataset object [BulkRNASeq Loaders Example](https://github.com/J-81/dp_tools/blob/main/dp_tools/bulkRNASeq/loaders.py)
    - Note: It is recommended to use the general data asset find function [source](dp_tools/bulkRNASeq/locaters.py#L20) coupled with a data asset yaml config [example](dp_tools/config/bulkRNASeq_v0.yaml)
1. Write checks [BulkRNASeq Checks Example](https://github.com/J-81/dp_tools/blob/main/dp_tools/bulkRNASeq/checks.py)
1. Incorporate those checks into a VV_Protocol [BulkRNASeq Protocol Example](https://github.com/J-81/dp_tools/blob/main/dp_tools/bulkRNASeq/vv_protocols.py)

# Additional these tests show how a top level VV script would likely look
1. First loading data [Excerpt Loading Test](https://github.com/J-81/dp_tools/blob/988a9d0a139404c61384479e035f307c4a23ceae/tests/test_loaders.py#L40)
1. Then running the protocol on the data (ignore the caplog stuff, that is part testing) [Excerpt Validation Test](https://github.com/J-81/dp_tools/blob/988a9d0a139404c61384479e035f307c4a23ceae/tests/test_validation.py#L86-L89)
