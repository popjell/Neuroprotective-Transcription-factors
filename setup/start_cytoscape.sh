#!/usr/bin/env bash
set -euo pipefail
export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
cytoscape_home="$HOME/lab/tools/cytoscape-unix-3.10.5"
if [[ ! -x "$cytoscape_home/cytoscape.sh" ]]; then
  printf 'Install Cytoscape with python3 setup/install_cytoscape.py first.\n' >&2
  exit 1
fi
exec "$cytoscape_home/cytoscape.sh" server "$@"
