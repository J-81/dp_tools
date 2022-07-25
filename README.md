# DP_TOOLS

## Installation

This software can be installed by cloning the repository and running:

``` bash
pip install </path/to/dp_tools>
```

Alternatively, released versions of this tool can be run in containers sourced from the associated [quay.io repository](https://quay.io/repository/j_81/dp_tools).

## Testing Setup

1. Install dp_tools (this codebase)
2. Install samtools (binary dependency for certain checks)
3. Install pytest (test runner)
4. Retrieve test data as follows:

    ``` bash
    git clone -b  NF_RCP-F  --depth=1 https://github.com/J-81/test-datasets-extended.git --single-branch && rm -rf test-datasets-extended/.git
    ```

5. export test data location variables

    ``` bash
    export TEST_ASSETS_DIR="/path/to/test-datasets-extended/testdata"
    ```

6. Run pytest againsts tests

    ``` bash
    pytest </path/to/dp_tools>
    ```

