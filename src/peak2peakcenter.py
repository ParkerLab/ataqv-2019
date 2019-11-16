import sys

peak = sys.argv[1]
peak = peak.replace(':', ' ').replace('-', ' ')
chrom, start, end = peak.split(' ')
start = int(start)
end = int(end)
center = int((end + start)/2)
center_end = center + 1
print(f'{chrom}\t{center}\t{center_end}')
