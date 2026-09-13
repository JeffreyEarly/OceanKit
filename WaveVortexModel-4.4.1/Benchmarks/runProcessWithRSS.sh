#!/bin/sh

sample_file="$1"
phase_file="$2"
interval="$3"
stdout_file="$4"
stderr_file="$5"
sampled_phase_file="${phase_file}.sampled"
shift 5

if [ "$1" != "--" ]; then
    printf '%s\n' "runProcessWithRSS.sh requires -- before the command." >&2
    exit 2
fi
shift

: > "$sample_file"
: > "$stdout_file"
: > "$stderr_file"
: > "$sampled_phase_file"

export WV_RSS_PHASE_ACK="$sampled_phase_file"
"$@" >"$stdout_file" 2>"$stderr_file" &
root_pid=$!
sample_index=0

sample_tree() {
    phase=$(sed -n '1p' "$phase_file" 2>/dev/null)
    values=$(ps -axo pid=,ppid=,rss= | awk -v root="$root_pid" '
        { pid[NR]=$1; parent[NR]=$2; rss[NR]=$3 }
        END {
            included[root]=1
            changed=1
            while (changed) {
                changed=0
                for (i=1; i<=NR; ++i) {
                    if (!included[pid[i]] && included[parent[i]]) {
                        included[pid[i]]=1
                        changed=1
                    }
                }
            }
            total=0; count=0
            for (i=1; i<=NR; ++i) {
                if (included[pid[i]]) { total+=rss[i]; count+=1 }
            }
            if (count>0) printf "%.0f\t%d", total, count
        }')
    if [ -n "$values" ]; then
        printf '%s\t%s\t%s\n' "$sample_index" "$phase" "$values" >> "$sample_file"
        printf '%s\n' "$phase" > "$sampled_phase_file"
    fi
    sample_index=$((sample_index + 1))
}

while kill -0 "$root_pid" 2>/dev/null; do
    sample_tree
    sleep "$interval"
done
wait "$root_pid"
exit_status=$?
sample_tree
exit "$exit_status"
