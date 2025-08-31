import csv
import os
import argparse
import math

def split_csv(input_file: str, output_dir: str, parts: int, encoding: str = 'utf-8') -> None:
    try:
        with open(input_file, newline='', encoding=encoding) as f:
            reader = list(csv.reader(f))
    except UnicodeDecodeError as e:
        raise UnicodeError(
            f"Failed to decode {input_file!r} with encoding {encoding}. "
            "Try specifying the correct encoding using --encoding."
        ) from e
    header = reader[0]
    rows = reader[1:]
    total_rows = len(rows)
    if parts <= 0:
        raise ValueError('Number of parts must be positive')
    chunk_size = math.ceil(total_rows / parts)
    os.makedirs(output_dir, exist_ok=True)
    for i in range(parts):
        start = i * chunk_size
        end = start + chunk_size
        chunk_rows = rows[start:end]
        if not chunk_rows:
            break
        part_path = os.path.join(output_dir, f"part_{i+1}.csv")
        with open(part_path, 'w', newline='', encoding=encoding) as out:
            writer = csv.writer(out)
            writer.writerow(header)
            writer.writerows(chunk_rows)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split CSV into N parts with header.')
    parser.add_argument('input_file', help='Path to input CSV file')
    parser.add_argument('output_dir', help='Directory where parts will be saved')
    parser.add_argument('parts', type=int, help='Number of parts to create')
    parser.add_argument('--encoding', default='utf-8',
                        help='File encoding (default: utf-8)')
    args = parser.parse_args()
    split_csv(args.input_file, args.output_dir, args.parts, args.encoding)
