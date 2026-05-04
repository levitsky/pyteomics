#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: rewrite_default_ns.sh [-p PREFIX] [-o OUTPUT] INPUT

Rewrite XML so elements in the document's default namespace are explicitly prefixed
(default prefix: ns0), then validate the result.

Options:
  -p PREFIX   Namespace prefix to use (default: ns0)
  -o OUTPUT   Output file path (default: INPUT with suffix .explicit-ns.xml)
  -h          Show this help
EOF
}

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Error: required command not found: $cmd" >&2
        exit 1
    fi
}

prefix="ns0"
output=""

while getopts ":p:o:h" opt; do
    case "$opt" in
        p)
            prefix="$OPTARG"
            ;;
        o)
            output="$OPTARG"
            ;;
        h)
            usage
            exit 0
            ;;
        :)
            echo "Error: option -$OPTARG requires an argument." >&2
            usage >&2
            exit 2
            ;;
        \?)
            echo "Error: invalid option -$OPTARG" >&2
            usage >&2
            exit 2
            ;;
    esac
done

shift $((OPTIND - 1))

if [[ $# -ne 1 ]]; then
    usage >&2
    exit 2
fi

input="$1"

if [[ ! -f "$input" ]]; then
    echo "Error: input file does not exist: $input" >&2
    exit 1
fi

if [[ ! "$prefix" =~ ^[A-Za-z_][A-Za-z0-9_.-]*$ ]]; then
    echo "Error: invalid XML prefix: $prefix" >&2
    exit 1
fi

if [[ -z "$output" ]]; then
    output="${input%.*}.explicit-ns.xml"
fi

require_cmd xmllint
require_cmd xmlstarlet
require_cmd mktemp

default_ns_uri="$(xmllint --xpath 'namespace-uri(/*)' "$input" 2>/dev/null || true)"
if [[ -z "$default_ns_uri" ]]; then
    echo "Error: input root element does not have a default namespace." >&2
    exit 1
fi

tmp_xslt="$(mktemp)"
trap 'rm -f "$tmp_xslt"' EXIT

cat > "$tmp_xslt" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:${prefix}="${default_ns_uri}"
    exclude-result-prefixes="${prefix}">
  <xsl:output method="xml" indent="yes"/>
  <xsl:strip-space elements="*"/>

  <xsl:template match="@*|text()|comment()|processing-instruction()">
    <xsl:copy/>
  </xsl:template>

  <xsl:template match="@*">
    <xsl:choose>
      <xsl:when test="namespace-uri() = '${default_ns_uri}'">
        <xsl:attribute name="${prefix}:{local-name()}" namespace="${default_ns_uri}">
          <xsl:value-of select="."/>
        </xsl:attribute>
      </xsl:when>
      <xsl:otherwise>
        <xsl:copy/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template match="*">
    <xsl:choose>
      <xsl:when test="namespace-uri() = '${default_ns_uri}'">
        <xsl:element name="${prefix}:{local-name()}" namespace="${default_ns_uri}">
          <xsl:apply-templates select="@*|node()"/>
        </xsl:element>
      </xsl:when>
      <xsl:otherwise>
        <xsl:copy>
          <xsl:apply-templates select="@*|node()"/>
        </xsl:copy>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>
</xsl:stylesheet>
EOF

xmlstarlet tr "$tmp_xslt" "$input" > "$output"

# Validation 1: well-formed XML
xmllint --noout "$output"

# Validation 2: root is prefixed with requested prefix
root_name="$(xmllint --xpath 'name(/*)' "$output" 2>/dev/null)"
if [[ "$root_name" != ${prefix}:* ]]; then
    echo "Error: output root element is not explicitly prefixed with '${prefix}'." >&2
    exit 1
fi

# Validation 3: no default namespace declaration on root
default_decl_count="$(xmllint --xpath 'count(/*/namespace::*[name()=""])' "$output" 2>/dev/null)"
if [[ "$default_decl_count" != "0" ]]; then
    echo "Error: output still has a default namespace declaration on the root element." >&2
    exit 1
fi

echo "Success: wrote '$output' with explicit namespace prefix '$prefix'."
echo "Validated: XML is well-formed, root is prefixed, default namespace declaration removed from root."