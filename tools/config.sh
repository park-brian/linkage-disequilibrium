# Shell twin of tools/config.mjs: source it, then `cd "$LD_DATA"` for data-relative paths.
#   . "$(dirname "$0")/../../tools/config.sh"     (from pipeline/<stage>/ or verify/)
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export REPO
export LD_DATA=${LD_DATA:-$(cd "$REPO/.." && pwd)/data}
export LD_RELEASE=${LD_RELEASE:-$LD_DATA/release_v4}
export LD_PANELS=${LD_PANELS:-$(cd "$REPO/.." && pwd)/db/data/1000G_population_samples}
export LD_VCF=${LD_VCF:-$(cd "$REPO/.." && pwd)/db/1000G}
export LDLINK_URL=${LDLINK_URL:-http://127.0.0.1:8090}
export LDLINK_HOME=${LDLINK_HOME:-$(cd "$REPO/.." && pwd)/ldlink-local}
mkdir -p "$LD_DATA"
