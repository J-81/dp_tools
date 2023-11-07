set -eux

OUT_DIR=$1

cp assets/isa_config.yaml ${OUT_DIR}
cp assets/extraction_conf.yaml ${OUT_DIR}/extraction_settings.yaml

# cp tests/assets/test.yaml ${OUT_DIR} # Disabled as test.yaml which contains user specified columns is not fully implemented
