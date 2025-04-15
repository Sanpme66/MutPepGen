from typing import List, Tuple, Dict, Set
import logging
from models.utils import dTime
import pandas as pd
import re
import os

logger = logging.getLogger(__name__)
# Constants
INIT_VARIABLES: Dict[str, int] = {
    'MAIN_counter': 0,
    'FILE_ERROR_counter': 0,
    'FILE_SUCCESS_counter': 0,
    'TRANSCRIPT_FOUND_counter': 0,
    'TRANSCRIPT_NOT_FOUND_counter': 0,
    'SUBSTITUTION_FOUND_counter': 0,
    'SUBSTITUTION_NOT_FOUND_counter': 0,
    'SUBSTITUTION_SUCCESS_counter': 0,
    'SUBSTITUTION_ERROR_counter': 0,
    'POSITION_FOUND_counter': 0,    
    'POSITION_2ND_ATTEMPT_FOUND_counter': 0,
    'POSITION_3RD_ATTEMPT_FOUND_counter': 0,
    'POSITION_NOT_FOUND_counter': 0,
    'UNIPROTtoGRch38_NOT_FOUND_counter': 0,
    'MULTI_SEQ_FOUND_counter': 0,
    'MULTI_SEQ_POSITION_FOUND_counter': 0,
    'SUBSTITUTION_FOUND_2ND_ATTEMPT_counter': 0,
    'SUBSTITUTION_FOUND_3RD_ATTEMPT_counter': 0,
}

class Files_Manager:
    def __init__(self):
        self.file = None
    def check_file(self):
        return self.file
    #https://stackoverflow.com/questions/40745686/python-process-file-using-multiple-cores
    @staticmethod
    def file_block(fp, number_of_blocks, block):
        '''
        A generator that splits a file into blocks and iterates
        over the lines of one of the blocks.

        '''

        assert 0 <= block and block < number_of_blocks
        assert 0 < number_of_blocks

        fp.seek(0,2)
        file_size = fp.tell()

        ini = file_size * block / number_of_blocks
        end = file_size * (1 + block) / number_of_blocks

        if ini <= 0:
            fp.seek(0)
        else:
            fp.seek(ini-1)
            fp.readline()

        while fp.tell() < end:
            yield fp.readline()
            
            
    def _find_fileType(self):
        if len(self.file.split('.')) > 1:
            return self.file.split('.')[-1]
        else:
            raise ValueError("File type not found")
    
    def check_permission(self):
        return self.file.mode
    
    def check_file_name(self):
        return self.file.name
    
    def check_file_encoding(self):
        return self.file.encoding
    
    def check_file_type(self):
        return self.file.mode
            
    def open_file(self, file_path):
        self.file = open(file_path, 'r')
        return self.file

    def read_file(self):
        return self.file.read()

    def close_file(self):
        self.file.close()
        self.file = None

class logs:
    def __init__(self):
        self.cwd = os.getcwd()
        self.logs = None
        self.dTime = dTime()
        
    def _init_logsFile(self):
        file_path = os.path.join(self.cwd, f'{self.dTime}_proccess.logs')
        logging.basicConfig(filename=file_path, filemode='w', format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
        
    def _add_info(self, info):
        logging.info(info)
    
    def _add_error(self, error):
        logging.error(error)
        
    def _add_warning(self, warning):
        logging.warning(warning)
        
    def _add_debug(self, debug):
        logging.debug(debug)

class UniProtParser:
    """
    Parser for UniProt database format to extract Ensembl transcript IDs
    and map them to protein sequences.
    """
    
    def __init__(self, log_callback=None):
        """
        Initialize the parser
        
        Args:
            log_callback: Function to call for logging messages
        """
        self.sequences = {}
        self.log = log_callback if log_callback else print
    
    def parse_file(self, file_path):
        """
        Parse a UniProt format file and extract Ensembl IDs and sequences
        
        Args:
            file_path: Path to the UniProt format file
            
        Returns:
            Dictionary of ENST IDs to sequences
        """
        self.log(f"Parsing UniProt format file: {file_path}")
        
        try:
            # Determine delimiter based on file extension
            file_ext = os.path.splitext(file_path)[1].lower()
            if file_ext in ['.tsv', '.txt']:
                delimiter = '\t'
            else:  # .csv
                delimiter = ','
            
            # Read the file
            df = pd.read_csv(file_path, delimiter=delimiter, low_memory=False)
            
            # First, identify the relevant columns
            ensembl_col = self._find_ensembl_column(df)
            sequence_col = self._find_sequence_column(df)
            
            if not ensembl_col:
                # Look for a column containing Ensembl IDs in the data
                for col in df.columns:
                    sample = df[col].dropna().astype(str).head(10)
                    if any('ENST' in str(val) for val in sample):
                        ensembl_col = col
                        break
            
            if not ensembl_col:
                self.log("No column containing Ensembl IDs found")
                return self.sequences
            
            if not sequence_col:
                # Try to find the sequence column by checking for long string content
                for col in df.columns:
                    if col != ensembl_col:  # Skip the Ensembl column
                        sample = df[col].dropna().astype(str).head(5)
                        if any(len(str(val)) > 50 and self._is_likely_protein_sequence(str(val)) for val in sample):
                            sequence_col = col
                            break
            
            if not sequence_col:
                self.log("No sequence column found")
                return self.sequences
            
            self.log(f"Using column '{ensembl_col}' for Ensembl IDs and '{sequence_col}' for sequences")
            
            # Process each row
            count = 0
            for idx, row in df.iterrows():
                if pd.isna(row[ensembl_col]) or pd.isna(row[sequence_col]):
                    continue
                
                # Get the sequence
                sequence = str(row[sequence_col]).strip()
                
                # Parse Ensembl IDs from the Ensembl column
                ensembl_text = str(row[ensembl_col])
                
                # For UniProt format, Ensembl IDs are often in quotes and semicolon-separated
                # Example: "ENST00000436697.3; ENSP00000484893.1; ENSG00000225973.4.";"ENST00000567948.1; ...
                
                # Extract all ENST IDs
                enst_ids = self._extract_all_enst_ids(ensembl_text)
                
                if enst_ids:
                    for enst_id in enst_ids:
                        # Remove version number if present
                        if "." in enst_id:
                            enst_id = enst_id.split(".")[0]
                        self.sequences[enst_id] = sequence
                        count += 1
            
            self.log(f"Successfully mapped {count} ENST IDs to sequences from UniProt format")
            
        except Exception as e:
            self.log(f"Error parsing UniProt format file: {str(e)}")
            raise
        
        return self.sequences
    
    def _find_ensembl_column(self, df):
        """Find the column containing Ensembl IDs"""
        ensembl_keywords = ['ensembl', 'enst', 'transcript']
        
        for col in df.columns:
            col_lower = str(col).lower()
            if any(keyword in col_lower for keyword in ensembl_keywords):
                return col
        
        return None
    
    def _find_sequence_column(self, df):
        """Find the column containing protein sequences"""
        sequence_keywords = ['sequence', 'seq', 'protein']
        
        for col in df.columns:
            col_lower = str(col).lower()
            if any(keyword in col_lower for keyword in sequence_keywords):
                return col
        
        return None
    
    def _extract_all_enst_ids(self, text):
        """Extract all ENST IDs from text"""
        # Handle UniProt format with quotes and semicolons
        # First, find all quoted sections
        quoted_sections = re.findall(r'"([^"]*)"', text)
        
        # For each section, find ENST IDs
        enst_ids = []
        
        # Process quoted sections
        for section in quoted_sections:
            matches = re.findall(r'(ENST\d+(?:\.\d+)?)', section)
            enst_ids.extend(matches)
        
        # Also look for ENSTs not in quotes
        matches = re.findall(r'(ENST\d+(?:\.\d+)?)', text)
        enst_ids.extend(matches)
        
        # Remove duplicates
        return list(set(enst_ids))
    
    def _is_likely_protein_sequence(self, text):
        """Check if a string is likely to be a protein sequence"""
        # Protein sequences typically contain only amino acid letters
        amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
        
        # Clean the string and convert to uppercase
        text = re.sub(r'\s', '', text).upper()
        
        # Check if at least 80% of characters are amino acids
        if len(text) == 0:
            return False
            
        amino_acid_count = sum(1 for char in text if char in amino_acids)
        return amino_acid_count / len(text) >= 0.8
    
    def save_to_fasta(self, output_path):
        """Save the parsed sequences to a FASTA file"""
        with open(output_path, 'w') as f:
            for enst_id, sequence in self.sequences.items():
                f.write(f">{enst_id}\n{sequence}\n")
                
        self.log(f"Saved {len(self.sequences)} sequences to {output_path}")


if __name__ == "__main__":
    parser = UniProtParser()
    sequences = parser.parse_file("uniprot_data.tsv")
    parser.save_to_fasta("uniprot_sequences.fasta")