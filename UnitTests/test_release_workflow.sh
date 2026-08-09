#!/usr/bin/env bash

set -euo pipefail

repository_root="$(cd "$(dirname "$0")/.." && pwd)"
temporary_root="$(mktemp -d)"
trap 'rm -rf "$temporary_root"' EXIT

initialize_repository_pair() {
    local working="$1"
    local remote="$2"
    git init --bare "$remote" >/dev/null
    git init -b main "$working" >/dev/null
    git -C "$working" config user.name Test
    git -C "$working" config user.email test@example.com
    git -C "$working" remote add origin "$remote"
}

commit_initial_package() {
    local working="$1"
    mkdir -p "$working/resources"
    printf '{"name":"FixturePackage","version":"1.2.3"}\n' > "$working/resources/mpackage.json"
    printf 'initial\n' > "$working/tracked.txt"
    git -C "$working" add .
    git -C "$working" commit -m initial >/dev/null
    git -C "$working" push -u origin main >/dev/null
}

source_root="$temporary_root/source"
source_remote="$temporary_root/source.git"
initialize_repository_pair "$source_root" "$source_remote"
commit_initial_package "$source_root"
ocean_root="$temporary_root/ocean"
ocean_remote="$temporary_root/ocean.git"
initialize_repository_pair "$ocean_root" "$ocean_remote"
printf 'OceanKit\n' > "$ocean_root/README.md"
git -C "$ocean_root" add README.md
git -C "$ocean_root" commit -m initial >/dev/null
git -C "$ocean_root" push -u origin main >/dev/null

preflight_path="$temporary_root/preflight.json"
SOURCE_REPOSITORY_ROOT="$source_root" \
BUMP_TYPE=none \
SHOULD_BUILD_DOCUMENTATION=false \
SHOULD_PACKAGE_FOR_DISTRIBUTION=false \
SHOULD_CHECK_DOCUMENTATION=false \
SHOULD_PROMOTE_UNRELEASED=false \
DOCUMENTATION_PACKAGE_SPECIFIER= \
DOCUMENTATION_CHECK_TASK=docs:check \
WORKFLOW_REF=workflow@ref \
WORKFLOW_SHA=abc123 \
PREFLIGHT_PATH="$preflight_path" \
    "$repository_root/tools/preflight_mpm_release.sh"
jq -e '.mode == "legacy" and .version == "1.2.3"' "$preflight_path" >/dev/null

non_distribution_metadata="$temporary_root/non-distribution-metadata.json"
jq -n \
    --arg sourceStartSha "$(git -C "$source_root" rev-parse HEAD)" \
    '{schemaVersion:1,mode:"legacy",packageName:"FixturePackage",oldVersion:"1.2.3",version:"1.2.3",bumpType:"none",sourceStartSha:$sourceStartSha,distributionCreated:false,snapshotFolder:null,exportPath:null,releaseBodyPath:null,authoringPaths:[]}' \
    > "$non_distribution_metadata"
source_remote_before="$(git --git-dir="$source_remote" rev-parse main)"
ocean_remote_before="$(git --git-dir="$ocean_remote" rev-parse main)"
"$repository_root/tools/publish_mpm_release.sh" "$non_distribution_metadata" "$source_root"
[[ "$(git --git-dir="$source_remote" rev-parse main)" == "$source_remote_before" ]]
[[ "$(git --git-dir="$ocean_remote" rev-parse main)" == "$ocean_remote_before" ]]

git -C "$source_root" tag v1.2.4
if SOURCE_REPOSITORY_ROOT="$source_root" \
    BUMP_TYPE=patch \
    SHOULD_BUILD_DOCUMENTATION=false \
    SHOULD_PACKAGE_FOR_DISTRIBUTION=false \
    SHOULD_CHECK_DOCUMENTATION=false \
    SHOULD_PROMOTE_UNRELEASED=false \
    DOCUMENTATION_PACKAGE_SPECIFIER= \
    DOCUMENTATION_CHECK_TASK=docs:check \
    WORKFLOW_REF=workflow@ref \
    WORKFLOW_SHA=abc123 \
    PREFLIGHT_PATH="$preflight_path" \
    "$repository_root/tools/preflight_mpm_release.sh" >/dev/null 2>&1; then
    printf 'Existing tag preflight unexpectedly succeeded.\n' >&2
    exit 1
fi
git -C "$source_root" tag -d v1.2.4 >/dev/null

mkdir "$ocean_root/FixturePackage-1.2.4"
if SOURCE_REPOSITORY_ROOT="$source_root" \
    OCEANKIT_PUBLISH_ROOT="$ocean_root" \
    BUMP_TYPE=patch \
    SHOULD_BUILD_DOCUMENTATION=false \
    SHOULD_PACKAGE_FOR_DISTRIBUTION=true \
    SHOULD_CHECK_DOCUMENTATION=false \
    SHOULD_PROMOTE_UNRELEASED=false \
    DOCUMENTATION_PACKAGE_SPECIFIER=ClassDocumentation@1.3.0 \
    DOCUMENTATION_CHECK_TASK=docs:check \
    WORKFLOW_REF=workflow@ref \
    WORKFLOW_SHA=abc123 \
    PREFLIGHT_PATH="$preflight_path" \
    "$repository_root/tools/preflight_mpm_release.sh" >/dev/null 2>&1; then
    printf 'Existing snapshot preflight unexpectedly succeeded.\n' >&2
    exit 1
fi
rmdir "$ocean_root/FixturePackage-1.2.4"

printf 'changed\n' > "$source_root/tracked.txt"
export_root="$temporary_root/export/FixturePackage-1.2.4"
mkdir -p "$export_root"
printf 'payload\n' > "$export_root/package.m"
release_body="$temporary_root/release-body.md"
printf '%s\n' '- Quotes: `code`, "double", 100%.' '' '- Second line.' > "$release_body"
metadata="$temporary_root/metadata.json"
jq -n \
    --arg sourceStartSha "$(git -C "$source_root" rev-parse HEAD)" \
    --arg exportPath "$export_root" \
    --arg releaseBodyPath "$release_body" \
    '{schemaVersion:1,mode:"pilot",packageName:"FixturePackage",oldVersion:"1.2.3",version:"1.2.4",bumpType:"patch",sourceStartSha:$sourceStartSha,distributionCreated:true,snapshotFolder:"FixturePackage-1.2.4",exportPath:$exportPath,releaseBodyPath:$releaseBodyPath,authoringPaths:["tracked.txt"]}' \
    > "$metadata"

stub_bin="$temporary_root/bin"
mkdir "$stub_bin"
cat > "$stub_bin/gh" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$@" > "$GH_ARGUMENTS"
while (( $# > 0 )); do
    if [[ "$1" == "--notes-file" ]]; then
        cp "$2" "$GH_NOTES_COPY"
        break
    fi
    shift
done
STUB
chmod +x "$stub_bin/gh"
trace="$temporary_root/trace"
gh_arguments="$temporary_root/gh-arguments"
gh_notes_copy="$temporary_root/gh-notes.md"
PATH="$stub_bin:$PATH" \
GH_ARGUMENTS="$gh_arguments" \
GH_NOTES_COPY="$gh_notes_copy" \
GITHUB_REF_NAME=main \
GITHUB_REPOSITORY=example/fixture \
OCEANKIT_RELEASE_TRACE="$trace" \
    "$repository_root/tools/publish_mpm_release.sh" "$metadata" "$source_root" "$ocean_root"

expected_trace=$'authoring-push\noceankit-push\ntag-push\ngithub-release'
[[ "$(cat "$trace")" == "$expected_trace" ]]
git --git-dir="$source_remote" show-ref --verify --quiet refs/tags/v1.2.4
git --git-dir="$ocean_remote" show main:FixturePackage-1.2.4/package.m >/dev/null
grep -Fx -- '--notes-file' "$gh_arguments" >/dev/null
grep -Fx -- "$release_body" "$gh_arguments" >/dev/null
cmp "$release_body" "$gh_notes_copy"

failed_source_root="$temporary_root/failed-source"
failed_source_remote="$temporary_root/failed-source.git"
initialize_repository_pair "$failed_source_root" "$failed_source_remote"
commit_initial_package "$failed_source_root"
failed_ocean_root="$temporary_root/failed-ocean"
failed_ocean_remote="$temporary_root/failed-ocean.git"
initialize_repository_pair "$failed_ocean_root" "$failed_ocean_remote"
printf 'OceanKit\n' > "$failed_ocean_root/README.md"
git -C "$failed_ocean_root" add README.md
git -C "$failed_ocean_root" commit -m initial >/dev/null
git -C "$failed_ocean_root" push -u origin main >/dev/null
cat > "$failed_ocean_remote/hooks/pre-receive" <<'HOOK'
#!/usr/bin/env bash
exit 1
HOOK
chmod +x "$failed_ocean_remote/hooks/pre-receive"
printf 'changed\n' > "$failed_source_root/tracked.txt"
failed_export="$temporary_root/failed-export/FixturePackage-1.2.4"
mkdir -p "$failed_export"
printf 'payload\n' > "$failed_export/package.m"
failed_metadata="$temporary_root/failed-metadata.json"
jq -n \
    --arg sourceStartSha "$(git -C "$failed_source_root" rev-parse HEAD)" \
    --arg exportPath "$failed_export" \
    --arg releaseBodyPath "$release_body" \
    '{schemaVersion:1,mode:"pilot",packageName:"FixturePackage",oldVersion:"1.2.3",version:"1.2.4",bumpType:"patch",sourceStartSha:$sourceStartSha,distributionCreated:true,snapshotFolder:"FixturePackage-1.2.4",exportPath:$exportPath,releaseBodyPath:$releaseBodyPath,authoringPaths:["tracked.txt"]}' \
    > "$failed_metadata"
failed_trace="$temporary_root/failed-trace"
if PATH="$stub_bin:$PATH" \
    GH_ARGUMENTS="$temporary_root/failed-gh-arguments" \
    GITHUB_REF_NAME=main \
    GITHUB_REPOSITORY=example/fixture \
    OCEANKIT_RELEASE_TRACE="$failed_trace" \
    "$repository_root/tools/publish_mpm_release.sh" "$failed_metadata" "$failed_source_root" "$failed_ocean_root" >/dev/null 2>&1; then
    printf 'Rejected OceanKit push unexpectedly succeeded.\n' >&2
    exit 1
fi
git --git-dir="$failed_source_remote" log -1 --format=%s | grep -Fx 'Bump version to v1.2.4' >/dev/null
if git --git-dir="$failed_source_remote" show-ref --verify --quiet refs/tags/v1.2.4; then
    printf 'Tag was created after a rejected OceanKit push.\n' >&2
    exit 1
fi
[[ ! -e "$temporary_root/failed-gh-arguments" ]]

[[ -z "$(git -C "$source_root" status --porcelain=v1 --untracked-files=all)" ]]
[[ -z "$(git -C "$ocean_root" status --porcelain=v1 --untracked-files=all)" ]]
[[ -z "$(git -C "$failed_source_root" status --porcelain=v1 --untracked-files=all)" ]]
[[ -z "$(git -C "$failed_ocean_root" status --porcelain=v1 --untracked-files=all)" ]]

workflow="$repository_root/.github/workflows/reusable-mpm-release.yml"
grep -F 'uses: actions/checkout@v6' "$workflow" >/dev/null
grep -F 'uses: matlab-actions/setup-matlab@v3' "$workflow" >/dev/null
grep -F 'uses: matlab-actions/run-command@v3' "$workflow" >/dev/null
grep -F 'repository: ${{ job.workflow_repository }}' "$workflow" >/dev/null
grep -F 'ref: ${{ job.workflow_sha }}' "$workflow" >/dev/null
grep -F 'if: ${{ inputs.shouldPackageForDistribution == '\''true'\'' }}' "$workflow" >/dev/null
if grep -F 'command:' "$workflow" | grep -F '${{ inputs.' >/dev/null; then
    printf 'A free-form workflow input is interpolated into MATLAB source.\n' >&2
    exit 1
fi

printf 'Release workflow shell fixtures passed.\n'
