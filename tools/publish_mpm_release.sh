#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    printf 'Usage: %s METADATA SOURCE_REPOSITORY [OCEANKIT_REPOSITORY]\n' "$0" >&2
    exit 2
fi
metadata_path="$1"
source_root="$2"
oceankit_root="${3:-}"

if [[ "$(jq -er '.schemaVersion' "$metadata_path")" != "1" ]]; then
    printf 'Unsupported release metadata schema.\n' >&2
    exit 1
fi
if ! jq -e '.authoringPaths | type == "array" and all(.[]; type == "string")' "$metadata_path" >/dev/null; then
    printf 'Metadata authoringPaths must be an array of relative paths.\n' >&2
    exit 1
fi

trace_event() {
    if [[ -n "${OCEANKIT_RELEASE_TRACE:-}" ]]; then
        printf '%s\n' "$1" >> "$OCEANKIT_RELEASE_TRACE"
    fi
}

mode="$(jq -er '.mode' "$metadata_path")"
version="$(jq -er '.version' "$metadata_path")"
bump_type="$(jq -er '.bumpType' "$metadata_path")"
distribution_created="$(jq -r '.distributionCreated' "$metadata_path")"
if [[ "$distribution_created" != "true" && "$distribution_created" != "false" ]]; then
    printf 'Metadata distributionCreated must be a boolean.\n' >&2
    exit 1
fi
snapshot_folder="$(jq -r '.snapshotFolder // empty' "$metadata_path")"
export_path="$(jq -r '.exportPath // empty' "$metadata_path")"
release_body_path="$(jq -r '.releaseBodyPath // empty' "$metadata_path")"
source_start_sha="$(jq -er '.sourceStartSha' "$metadata_path")"
package_name="$(jq -er '.packageName' "$metadata_path")"

case "$mode" in
    legacy|pilot) ;;
    *) printf 'Unsupported release mode: %s\n' "$mode" >&2; exit 1 ;;
esac
case "$bump_type" in
    none|patch|minor|major) ;;
    *) printf 'Unsupported bump type: %s\n' "$bump_type" >&2; exit 1 ;;
esac

if [[ "$(git -C "$source_root" rev-parse HEAD)" != "$source_start_sha" ]]; then
    printf 'Source checkout moved after preflight; refusing to publish.\n' >&2
    exit 1
fi

git -C "$source_root" config user.name "github-actions[bot]"
git -C "$source_root" config user.email "github-actions[bot]@users.noreply.github.com"
authoring_paths=()
while IFS= read -r authoring_path; do
    case "$authoring_path" in
        ""|/*|../*|*/../*|*/..) printf 'Unsafe authoring path: %s\n' "$authoring_path" >&2; exit 1 ;;
    esac
    authoring_paths+=("$authoring_path")
done < <(jq -r '.authoringPaths[]?' "$metadata_path") || true
if (( ${#authoring_paths[@]} > 0 )); then
    git -C "$source_root" add -- "${authoring_paths[@]}"
fi
if ! git -C "$source_root" diff --cached --quiet; then
    git -C "$source_root" commit -m "Bump version to v$version"
    source_branch="${GITHUB_REF_NAME:-$(git -C "$source_root" branch --show-current)}"
    trace_event authoring-push
    git -C "$source_root" push origin "HEAD:$source_branch"
fi

if [[ "$distribution_created" == "true" ]]; then
    if [[ -z "$oceankit_root" || ! -d "$oceankit_root/.git" ]]; then
        printf 'Distribution metadata requires a writable OceanKit checkout.\n' >&2
        exit 1
    fi
    if [[ -z "$snapshot_folder" || ! -d "$export_path" ]]; then
        printf 'Distribution metadata does not identify a valid export.\n' >&2
        exit 1
    fi
    target_path="$oceankit_root/$snapshot_folder"
    git -C "$oceankit_root" fetch origin main
    if git -C "$oceankit_root" cat-file -e "origin/main:$snapshot_folder" 2>/dev/null; then
        printf 'OceanKit snapshot already exists on origin/main: %s\n' "$snapshot_folder" >&2
        exit 1
    fi
    if [[ "$mode" == "pilot" && -e "$target_path" ]]; then
        printf 'OceanKit snapshot already exists: %s\n' "$snapshot_folder" >&2
        exit 1
    fi
    mkdir -p "$target_path"
    rsync -a "$export_path/" "$target_path/"
    git -C "$oceankit_root" config user.name "github-actions[bot]"
    git -C "$oceankit_root" config user.email "github-actions[bot]@users.noreply.github.com"
    git -C "$oceankit_root" add -- "$snapshot_folder"
    if git -C "$oceankit_root" diff --cached --quiet; then
        printf 'No OceanKit snapshot changes were produced.\n' >&2
        exit 1
    fi
    git -C "$oceankit_root" commit -m "Add $snapshot_folder (from ${GITHUB_REPOSITORY:-unknown}@$source_start_sha)"
    trace_event oceankit-push
    git -C "$oceankit_root" push origin HEAD:main
fi

if [[ "$bump_type" != "none" ]]; then
    tag="v$version"
    if git -C "$source_root" show-ref --verify --quiet "refs/tags/$tag" || \
            git -C "$source_root" ls-remote --exit-code --tags origin "refs/tags/$tag" >/dev/null 2>&1; then
        printf 'Release tag already exists: %s\n' "$tag" >&2
        exit 1
    fi
    trace_event tag-push
    git -C "$source_root" tag "$tag"
    git -C "$source_root" push origin "$tag"
    trace_event github-release
    if [[ -n "$release_body_path" ]]; then
        gh release create "$tag" --repo "${GITHUB_REPOSITORY:?}" --title "$tag" --notes-file "$release_body_path"
    else
        gh release create "$tag" --repo "${GITHUB_REPOSITORY:?}" --title "$tag" --generate-notes
    fi
fi

printf 'Published release candidate for %s %s.\n' "$package_name" "$version"
