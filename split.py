import sys

def split_fasta(fasta_file, output_prefix, chunk_size=500*(10**6)):
    chunks = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if 'current_chunk' in locals():
                    chunks.append((header, ''.join(current_chunk)))
                header = line.strip()
                current_chunk = []
            else:
                current_chunk.append(line.strip())
        if 'current_chunk' in locals():  # Add the last chromosome chunk
            chunks.append((header, ''.join(current_chunk)))

    sub_chunks = {}
    for header, sequence in chunks:
        chr_id = header[1:]
        len_seq = len(sequence)
        num_sub_chunks = max(1, len_seq // chunk_size + (1 if len_seq % chunk_size else 0))
        sub_chunk_size = int(len_seq / num_sub_chunks)
        for i in range(num_sub_chunks):
            sub_chunk_id = '{}_{}'.format(chr_id, i + 1)
            start = i * sub_chunk_size
            end = start + sub_chunk_size if i < num_sub_chunks - 1 else len_seq
            sub_chunks[sub_chunk_id] = (header, start + 1, end)  # convert to 1-based index
            with open('{}_{}.fasta'.format(output_prefix, sub_chunk_id), 'w') as out:
                out.write('>{}\n'.format(sub_chunk_id))
                out.write(sequence[start:end] + '\n')
    return sub_chunks

def update_gff(gff_file, sub_chunks, output_file):
    with open(gff_file, 'r') as gff, open(output_file, 'w') as out:
        for line in gff:
            if line.startswith("#"):
                out.write(line)
                continue
            parts = line.strip().split('\t')
            chr_id = parts[0]
            start, end = int(parts[3]), int(parts[4])
            for sub_chunk_id in sub_chunks:
                _header, sub_start, sub_end = sub_chunks[sub_chunk_id]
                if chr_id in _header:
                    if sub_start <= start <= sub_end:
                        parts[0] = sub_chunk_id
                        parts[3] = str(start - sub_start + 1)  # convert to 0-based index
                        parts[4] = str(end - sub_start + 1)  # keep within the new sub-chunk
                        break
            out.write('\t'.join(parts) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python2.7 split.py <fasta_file> <gff_file> <output_prefix> <gff_output_file>")
        sys.exit()

    fasta_file = sys.argv[1]
    gff_file = sys.argv[2]
    output_prefix = sys.argv[3]
    gff_output_file = sys.argv[4]

    sub_chunks = split_fasta(fasta_file, output_prefix)
    update_gff(gff_file, sub_chunks, gff_output_file)
