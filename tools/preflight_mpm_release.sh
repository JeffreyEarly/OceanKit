#!/usr/bin/env bash

set -euo pipefail

required_variables=(
    SOURCE_REPOSITORY_ROOT
    BUMP_TYPE
    SHOULD_BUILD_DOCUMENTATION
    SHOULD_PACKAGE_FOR_DISTRIBUTION
    SHOULD_CHECK_DOCUMENTATION
    SHOULD_PROMOTE_UNRELEASED
    DOCUMENTATION_PACKAGE_SPECIFIER
    DOCUMENTATION_CHECK_TASK
    WORKFLOW_REF
    WORKFLOW_SHA
    PREFLIGHT_PATH
)
for variable_name in "${required_variables[@]}"; do
    if [[ -z "${!variable_name+x}" ]]; then
        printf 'Missing required preflight variable: %s\n' "$variable_name" >&2
        exit 1
    fi
done

case "$BUMP_TYPE" in
    none|patch|minor|major) ;;
    *) printf 'Invalid bump type: %s\n' "$BUMP_TYPE" >&2; exit 1 ;;
esac
for boolean_name in SHOULD_BUILD_DOCUMENTATION SHOULD_PACKAGE_FOR_DISTRIBUTION SHOULD_CHECK_DOCUMENTATION SHOULD_PROMOTE_UNRELEASED; do
    boolean_value="${!boolean_name}"
    if [[ "$boolean_value" != "true" && "$boolean_value" != "false" ]]; then
        printf '%s must be true or false, not %s.\n' "$boolean_name" "$boolean_value" >&2
        exit 1
    fi
done

manifest_path="$SOURCE_REPOSITORY_ROOT/resources/mpackage.json"
if [[ ! -f "$manifest_path" ]]; then
    printf 'Package manifest not found: %s\n' "$manifest_path" >&2
    exit 1
fi
package_name="$(jq -er '.name | select(type == "string" and length > 0)' "$manifest_path")"
old_version="$(jq -er '.version | select(type == "string" and test("^[0-9]+\\.[0-9]+\\.[0-9]+$"))' "$manifest_path")"
IFS=. read -r major minor patch <<< "$old_version"
case "$BUMP_TYPE" in
    major) version="$((major + 1)).0.0" ;;
    minor) version="$major.$((minor + 1)).0" ;;
    patch) version="$major.$minor.$((patch + 1))" ;;
    none) version="$old_version" ;;
esac

pilot_enabled=false
if [[ -n "$DOCUMENTATION_PACKAGE_SPECIFIER" || "$SHOULD_CHECK_DOCUMENTATION" == "true" || "$SHOULD_PROMOTE_UNRELEASED" == "true" ]]; then
    pilot_enabled=true
fi

if [[ -n "$(git -C "$SOURCE_REPOSITORY_ROOT" status --porcelain=v1 --untracked-files=all)" ]]; then
    printf 'Release preflight requires a clean source checkout.\n' >&2
    git -C "$SOURCE_REPOSITORY_ROOT" status --short >&2
    exit 1
fi

if [[ "$BUMP_TYPE" != "none" ]] && git -C "$SOURCE_REPOSITORY_ROOT" show-ref --verify --quiet "refs/tags/v$version"; then
    printf 'Release tag already exists: v%s\n' "$version" >&2
    exit 1
fi

if [[ "$pilot_enabled" == "true" && "$SHOULD_CHECK_DOCUMENTATION" == "true" && -z "$DOCUMENTATION_PACKAGE_SPECIFIER" ]]; then
    printf 'Documentation checking requires an exact documentation package specifier.\n' >&2
    exit 1
fi
if [[ "$pilot_enabled" == "true" && "$SHOULD_CHECK_DOCUMENTATION" == "true" ]]; then
    if [[ -z "$DOCUMENTATION_CHECK_TASK" || ! -f "$SOURCE_REPOSITORY_ROOT/buildfile.m" ]]; then
        printf 'Documentation checking requires a task name and buildfile.m.\n' >&2
        exit 1
    fi
fi
if [[ "$pilot_enabled" == "true" && "$SHOULD_BUILD_DOCUMENTATION" == "true" && ! -f "$SOURCE_REPOSITORY_ROOT/tools/build_website_documentation.m" ]]; then
    printf 'Requested documentation builder is missing.\n' >&2
    exit 1
fi
if [[ "$SHOULD_PROMOTE_UNRELEASED" == "true" ]]; then
    if [[ "$BUMP_TYPE" == "none" || "$SHOULD_BUILD_DOCUMENTATION" != "true" ]]; then
        printf 'Unreleased promotion requires a version bump and documentation build.\n' >&2
        exit 1
    fi
    if [[ ! -f "$SOURCE_REPOSITORY_ROOT/CHANGELOG.md" ]]; then
        printf 'Unreleased promotion requires CHANGELOG.md.\n' >&2
        exit 1
    fi
fi

snapshot_folder="$package_name-$version"
if [[ "$pilot_enabled" == "true" && "$SHOULD_PACKAGE_FOR_DISTRIBUTION" == "true" ]]; then
    if [[ -z "${OCEANKIT_PUBLISH_ROOT:-}" || ! -d "$OCEANKIT_PUBLISH_ROOT/.git" ]]; then
        printf 'Distribution preflight requires a writable OceanKit checkout.\n' >&2
        exit 1
    fi
    if [[ -e "$OCEANKIT_PUBLISH_ROOT/$snapshot_folder" ]]; then
        printf 'OceanKit snapshot already exists: %s\n' "$snapshot_folder" >&2
        exit 1
    fi
fi

release_date="$(date -u +%F)"
mkdir -p "$(dirname "$PREFLIGHT_PATH")"
temporary_path="$PREFLIGHT_PATH.tmp"
jq -n \
    --argjson schemaVersion 1 \
    --arg mode "$([[ "$pilot_enabled" == "true" ]] && printf pilot || printf legacy)" \
    --arg workflowRef "$WORKFLOW_REF" \
    --arg workflowSha "$WORKFLOW_SHA" \
    --arg sourceStartSha "$(git -C "$SOURCE_REPOSITORY_ROOT" rev-parse HEAD)" \
    --arg packageName "$package_name" \
    --arg oldVersion "$old_version" \
    --arg version "$version" \
    --arg bumpType "$BUMP_TYPE" \
    --arg releaseDate "$release_date" \
    --arg snapshotFolder "$snapshot_folder" \
    '{schemaVersion:$schemaVersion,mode:$mode,workflowRef:$workflowRef,workflowSha:$workflowSha,sourceStartSha:$sourceStartSha,packageName:$packageName,oldVersion:$oldVersion,version:$version,bumpType:$bumpType,releaseDate:$releaseDate,snapshotFolder:$snapshotFolder}' \
    > "$temporary_path"
mv "$temporary_path" "$PREFLIGHT_PATH"

printf 'Release preflight: %s %s -> %s (%s)\n' "$package_name" "$old_version" "$version" "$([[ "$pilot_enabled" == "true" ]] && printf pilot || printf legacy)"
if [[ -n "${GITHUB_STEP_SUMMARY:-}" ]]; then
    {
        printf '### MPM release preflight\n\n'
        printf -- '- Workflow: `%s` (`%s`)\n' "$WORKFLOW_REF" "$WORKFLOW_SHA"
        printf -- '- Package: `%s`\n' "$package_name"
        printf -- '- Version: `%s` → `%s`\n' "$old_version" "$version"
        printf -- '- Mode: `%s`\n' "$([[ "$pilot_enabled" == "true" ]] && printf pilot || printf legacy)"
    } >> "$GITHUB_STEP_SUMMARY"
fi
