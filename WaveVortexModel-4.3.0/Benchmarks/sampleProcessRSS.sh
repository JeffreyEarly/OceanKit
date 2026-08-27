#!/bin/sh

pid="$1"
phase_file="$2"
stop_file="$3"
output_file="$4"
interval="$5"
sample_index=0

while [ ! -e "$stop_file" ]; do
    phase=$(sed -n '1p' "$phase_file" 2>/dev/null)
    rss_kib=$(ps -o rss= -p "$pid" 2>/dev/null | tr -d ' ')
    if [ -n "$rss_kib" ]; then
        printf '%s\t%s\t%s\n' "$sample_index" "$phase" "$rss_kib" >> "$output_file"
    fi
    sample_index=$((sample_index + 1))
    sleep "$interval"
done
