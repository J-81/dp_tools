# DP_TOOLS

## Installation

This software can be installed by cloning the repository and running:

``` bash
pip install </path/to/dp_tools>
```

Alternatively, released versions of this tool can be run in containers sourced from the associated [quay.io repository]().

## Testing Setup

1. Install dp_tools
2. Install samtools
3. Retrieve test data as follows:

    ``` bash
    git clone -b  NF_RCP-F  --depth=1 https://github.com/J-81/test-datasets-extended.git --single-branch && rm -rf test-datasets-extended/.git
    ```

4. export test data location variables

    ``` bash
    export TEST_ASSETS_DIR="/path/to/test-datasets-extended/testdata"
    ```
