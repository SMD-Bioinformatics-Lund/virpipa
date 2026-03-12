#!/usr/bin/env python3
"""
Standalone script to download the hbv_result table from geno2pheno HCV rules page.
Saves to properly formatted CSV with proper quoting for fields containing commas.
"""

import requests
from bs4 import BeautifulSoup
import csv
import sys
import os
from datetime import datetime

def escape_csv_field(field):
    """Escape a CSV field - csv.writer handles quoting automatically."""
    if field is None:
        return ''
    return str(field).strip()


def download_geno2pheno_rules(output_file=None, force_refresh=False):
    """
    Download the hbv_result table from geno2pheno rules page.

    Args:
        output_file: Path to output CSV file (default: hbv_result_rules.csv)
        force_refresh: Always download even if file exists
    """

    url = 'https://hcv.geno2pheno.org/index.php?page=Rules'
    default_output = 'hbv_result_rules.csv'

    # Set output file
    if output_file is None:
        output_file = default_output

    # Check if file exists and skip download if not force_refresh
    if not force_refresh and os.path.exists(output_file):
        print(f"CSV file already exists: {output_file}")
        print("Use --force-refresh to re-download")
        return output_file

    print(f"Fetching hbv_result table from {url}...")
    print(f"Output: {output_file}")
    print("-" * 60)

    try:
        # Fetch the page
        response = requests.get(url, timeout=60)
        response.raise_for_status()

        # Parse HTML
        soup = BeautifulSoup(response.content, 'html.parser')

        # Find the hbv_result table
        table = soup.find('table', {'class': 'hbv_result'})

        if not table:
            print("ERROR: Could not find hbv_result table")
            return None

        # Extract table data
        rows = []
        headers = []

        # Parse the header row - extract text from each cell individually
        # because the HTML has nested <th> tags with multiple <div> elements
        first_row = table.find('tr')
        if first_row:
            header_cells = first_row.find_all(['td', 'th'])
        else:
            header_cells = []

        headers = []
        for cell in header_cells:
            # Find all div tags and take the text from the first one (not the last)
            divs = cell.find_all('div')
            if divs:
                header_text = divs[0].get_text(strip=True)
            else:
                header_text = cell.get_text(strip=True)
            headers.append(header_text)

        # If headers are still empty, try manual mapping
        if not headers or len(headers) < 7:
            headers = ['Drugs', 'Region', 'Rule Definitions', 'Rule for Subtype',
                      'Drug licensed for genotype', 'Predictions', 'Reference']

        # Parse the actual data rows
        data_rows = table.find_all('tr')[1:]

        for row in data_rows:
            cells = row.find_all(['td', 'th'])
            row_data = []

            # Extract text from each cell individually
            for cell in cells:
                # Find all div tags and take the text from the first one (not the last)
                divs = cell.find_all('div')
                if divs:
                    cell_text = divs[0].get_text(strip=True)
                else:
                    cell_text = cell.get_text(strip=True)
                row_data.append(cell_text)

            # Only add if we got reasonable data (7 columns)
            if len(row_data) == 7:
                rows.append(row_data)

        # Save to CSV with proper quoting
        print(f"\nExtracted {len(rows)} rule entries")
        print(f"Saving to {output_file}...")

        drugs = set(r[0] for r in rows)
        regions = set(r[1] for r in rows)
        predictions = set(r[5] for r in rows)
        subtypes = set(r[3] for r in rows)

        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            # Write comments first
            csvfile.write(f"# Downloaded from: {url}\n")
            csvfile.write(f"# Last updated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            csvfile.write("#\n")

            # Write header row with proper escaping
            writer.writerow([escape_csv_field(h) for h in headers])

            # Write data rows with proper escaping
            for row in rows:
                writer.writerow([escape_csv_field(field) for field in row])

        print("-" * 60)
        print("SUMMARY:")
        print(f"  Total rule entries: {len(rows)}")
        print(f"  Unique drugs: {len(drugs)}")
        print(f"  Unique regions: {len(regions)}")
        print(f"  Unique predictions: {len(predictions)}")
        print(f"  Unique subtypes: {len(subtypes)}")
        print(f"  Drugs: {', '.join(sorted(drugs))}")
        print("-" * 60)
        print(f"\n✓ Saved successfully to {output_file}")
        print(f"\nColumns: {', '.join(headers)}")

        return output_file

    except requests.exceptions.RequestException as e:
        print(f"ERROR: Failed to download page: {e}")
        return None
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """Main function for command-line usage"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Download hbv_result table from geno2pheno HCV rules page',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                           # Download to default file (hbv_result_rules.csv)
  %(prog)s -o my_rules.csv           # Save to custom file
  %(prog)s --force-refresh           # Re-download even if file exists
  %(prog)s -o my_rules.csv --force-refresh
        """
    )

    parser.add_argument('-o', '--output', metavar='FILE',
                        help='Output CSV file (default: hbv_result_rules.csv)')
    parser.add_argument('--force-refresh', action='store_true',
                        help='Always download, even if file already exists')

    args = parser.parse_args()

    result = download_geno2pheno_rules(args.output, args.force_refresh)

    if result:
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()
