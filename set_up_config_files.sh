set -eux

OUT_DIR=$1

cp tests/assets/isa_config.yaml ${OUT_DIR}
cp tests/assets/config.yaml ${OUT_DIR}/extraction_settings.yaml

# cp tests/assets/test.yaml ${OUT_DIR} # Disabled as test.yaml which contains user specified columns is not fully implemented