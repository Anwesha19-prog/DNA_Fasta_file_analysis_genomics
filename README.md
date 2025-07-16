DNA Fasta File Analysis – Genomics

AUTHOR: Anwesha Sarkar
🧬 DNA Sequence Analysis
This repository contains a Python script for analyzing DNA sequences from FASTA files. The script identifies Open Reading Frames (ORFs) across reading frames and performs various sequence analyses, including repeat detection.

🚀 Features
ORF Identification: Detects the longest Open Reading Frames in specified reading frames.

Sequence Statistics: Calculates the longest and shortest sequence lengths.

Repeat Counting: Identifies the most frequent k-mers of specified lengths (e.g., 6, 7, 12).

🛠️ Installation
1. Clone the repository
bash
Copy
Edit
git clone https://github.com/Anwesha19-prog/dna-sequence-analysis.git
cd dna-sequence-analysis
2. Install dependencies
Make sure Python is installed, then install Biopython:

bash
Copy
Edit
pip install biopython
▶️ Usage
Prepare your FASTA file
Place your FASTA file in the directory or update the script with its path.

Run the script

Edit the script to include your FASTA file path:

python
Copy
Edit
fasta_file = r"your_fasta_file.fasta"  # replace with your file path
results = analyze_fasta(fasta_file)
Then run the script:

bash
Copy
Edit
python analyze_fasta.py
View Results
The script will output various sequence statistics. For example:

python
Copy
Edit
print(f"Number of records: {results['record_count']}")
print(f"Length of the longest sequence: {results['longest_seq_len']}")
print(f"Length of the shortest sequence: {results['shortest_seq_len']}")
print(f"Length of the longest ORF in reading frame 2: {results['longest_orf_rf2']}")
print(f"Starting position of the longest ORF in reading frame 3: {results['longest_orf_rf3_start']}")
print(f"Length of the longest ORF in reading frame 3: {results['longest_orf_rf3_len']}")
print(f"Length of the longest ORF in any forward reading frame: {results['longest_orf_any_frame']}")

📄 License
This project is licensed under the MIT License. See the LICENSE file for details.

MIT License

Copyright (c) 2025 Anwesha19-prog

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
