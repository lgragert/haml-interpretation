import xml.etree.ElementTree as ET
import csv
import os
import sys

# Define the namespace for HAML
HAML_NAMESPACE_URI = "urn:HAML.Namespace"
NAMESPACE = {'hamlns': HAML_NAMESPACE_URI}

# Register the namespace with a prefix (or as default) for output serialization
# Using '' as prefix makes it the default namespace in the output XML for elements from this URI
ET.register_namespace('', HAML_NAMESPACE_URI)

def create_namespaced_element(tag_name):
    """Helper function to create an element with the HAML namespace."""
    return ET.Element(f"{{{HAML_NAMESPACE_URI}}}{tag_name}")

def find_elements_ns(parent, path):
    """Helper function to find multiple elements using the HAML namespace."""
    return parent.findall(path, NAMESPACE)

def find_element_ns(parent, path):
    """Helper function to find a single element using the HAML namespace."""
    return parent.find(path, NAMESPACE)

def indent_xml(elem, level=0):
    """Helper function to add proper indentation to XML elements."""
    i = "\n" + level * "    "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "    "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for subelem in elem:
            indent_xml(subelem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def load_antigen_conversion_map(csv_file):
    """Load antigen conversion rules from CSV file."""
    conversion_map = {}
    try:
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Get two-field version of the IMGT_HLA_Allele
                allele = row['IMGT_HLA_Allele'].split(':')[:2]
                allele = ':'.join(allele)
                conversion_map[allele] = row['Antigen']
        return conversion_map
    except FileNotFoundError:
        print(f"Error: Antigen conversion file '{csv_file}' not found.")
        return None
    except Exception as e:
        print(f"Error reading antigen conversion file: {e}")
        return None

def convert_allele_to_antigen(allele, conversion_map):
    """Convert an HLA allele to its corresponding antigen."""
    # Handle heterodimer case
    if '&' in allele:
        parts = allele.split('&')
        # Only convert DQB1/DPB1 alleles in heterodimers
        for part in parts:
            if part.startswith('DQB1*') or part.startswith('DPB1*'):
                return conversion_map.get(part, part)
        return allele  # Return original if no DQB1/DPB1 found
    
    return conversion_map.get(allele, allele)  # Handle non-heterodimer case

def create_abl_string(pl_string, conversion_map):
    """Convert a PL string of alleles to an ABL string of antigens."""
    if not pl_string:
        return ""
    
    alleles = pl_string.split('+')
    antigens = []
    
    for allele in alleles:
        allele = allele.strip()
        antigen = convert_allele_to_antigen(allele, conversion_map)
        if antigen not in antigens:  # Avoid duplicates
            antigens.append(antigen)
            
    return '+'.join(antigens)

def get_bead_classification(mfi):
    """Helper function to determine bead classification based on MFI value."""
    if mfi > 2000:
        return "Positive", "8"
    elif 1000 <= mfi <= 2000:
        return "Borderline", "4"
    else:  # mfi < 1000
        return "Negative", "1"

def process_xml_file(input_file_path, output_file_path):
    """
    Parses the HAML XML file, adds interpretation blocks, and writes the modified XML.
    """
    # Load antigen conversion map
    conversion_map = load_antigen_conversion_map('antigen_conversion_table_with_rules.csv')
    if conversion_map is None:
        return

    try:
        tree = ET.parse(input_file_path)
        root = tree.getroot()
    except ET.ParseError as e:
        print(f"Error parsing XML file: {e}")
        return
    except FileNotFoundError:
        print(f"Error: Input file '{input_file_path}' not found.")
        return

    # Iterate through each patient
    for patient in find_elements_ns(root, 'hamlns:patient'):
        # Iterate through each sample within the patient
        for sample in find_elements_ns(patient, 'hamlns:sample'):
            # Iterate through each assay within the sample
            for assay in find_elements_ns(sample, 'hamlns:assay'):
                positive_hla_types = []
                questionable_hla_types = []
                negative_hla_types = []

                working_sample = find_element_ns(assay, 'hamlns:working-sample')
                if working_sample is None:
                    continue

                solid_phase_panel = find_element_ns(working_sample, 'hamlns:solid-phase-panel')
                if solid_phase_panel is None:
                    continue

                # Add interpretation software information to solid-phase-panel
                interpretation_software = create_namespaced_element('interpretation-software')
                interpretation_software.text = "Python-Demo"
                solid_phase_panel.append(interpretation_software)

                interpretation_software_version = create_namespaced_element('interpretation-software-version')
                interpretation_software_version.text = "0.0.1"
                solid_phase_panel.append(interpretation_software_version)

                # Initialize lists before bead processing
                positive_hla_types = []
                questionable_hla_types = []
                negative_hla_types = []

                # Process each bead
                for bead in find_elements_ns(solid_phase_panel, 'hamlns:bead'):
                    bead_info = find_element_ns(bead, 'hamlns:bead-info')
                    hla_target_type_el = find_element_ns(bead_info, 'hamlns:HLA-target-type')
                    
                    if hla_target_type_el is not None and hla_target_type_el.text is not None:
                        hla_target_type = hla_target_type_el.text
                        raw_data = find_element_ns(bead, 'hamlns:raw-data')
                        
                        # Get or create converted-data element
                        converted_data = find_element_ns(bead, 'hamlns:converted-data')
                        if converted_data is None:
                            converted_data = create_namespaced_element('converted-data')
                            bead.append(converted_data)

                        # Get MFI value and create Python-Demo interpretation
                        if raw_data is not None:
                            sample_raw_mfi_el = find_element_ns(raw_data, 'hamlns:sample-raw-MFI')
                            if sample_raw_mfi_el is not None and sample_raw_mfi_el.text is not None:
                                try:
                                    mfi = int(sample_raw_mfi_el.text)
                                    
                                    # Create new bead interpretation for Python-Demo
                                    bead_interpretation = create_namespaced_element('bead-interpretation')
                                    
                                    # Add classification entity
                                    classification_entity = create_namespaced_element('classification-entity')
                                    classification_entity.text = "Python-Demo"
                                    bead_interpretation.append(classification_entity)
                                    
                                    # Add interpretation reason
                                    interpretation_reason = create_namespaced_element('interpretation-reason')
                                    interpretation_reason.text = "MFI-cutoff-setting"
                                    bead_interpretation.append(interpretation_reason)
                                    
                                    # Add bead classification and rank based on MFI
                                    classification, rank = get_bead_classification(mfi)
                                    
                                    bead_classification = create_namespaced_element('bead-classification')
                                    bead_classification.text = classification
                                    bead_interpretation.append(bead_classification)
                                    
                                    bead_rank = create_namespaced_element('bead-rank')
                                    bead_rank.text = rank
                                    bead_interpretation.append(bead_rank)
                                    
                                    # Add to converted-data
                                    converted_data.append(bead_interpretation)
                                    
                                    # Add HLA type to appropriate list based on classification
                                    if classification == "Positive":
                                        positive_hla_types.append(hla_target_type)
                                    elif classification == "Borderline":
                                        questionable_hla_types.append(hla_target_type)
                                    elif classification == "Negative":
                                        negative_hla_types.append(hla_target_type)
                                        
                                except ValueError:
                                    print(f"Warning: Could not parse MFI value for bead")
                                    continue

                # Create the interpretation element
                interpretation_el = create_namespaced_element('interpretation')

                # Add interpretation software information
                interpretation_software = create_namespaced_element('interpretation-software')
                interpretation_software.text = "Python-Demo"
                interpretation_el.append(interpretation_software)

                interpretation_software_version = create_namespaced_element('interpretation-software-version')
                interpretation_software_version.text = "0.0.1"
                interpretation_el.append(interpretation_software_version)

                # Positive specificities
                positive_specificities_el = create_namespaced_element('positive-specificities')
                if positive_hla_types:
                    # Add PL string
                    hla_plstring_pos = create_namespaced_element('HLA-plstring')
                    pl_string = '+'.join(positive_hla_types)
                    hla_plstring_pos.text = pl_string
                    positive_specificities_el.append(hla_plstring_pos)
                    
                    # Add ABL string
                    hla_ablstring_pos = create_namespaced_element('HLA-ablstring')
                    hla_ablstring_pos.text = create_abl_string(pl_string, conversion_map)
                    positive_specificities_el.append(hla_ablstring_pos)
                interpretation_el.append(positive_specificities_el)

                # Questionable specificities
                questionable_specificities_el = create_namespaced_element('questionable-specificities')
                if questionable_hla_types:
                    # Add PL string
                    hla_plstring_ques = create_namespaced_element('HLA-plstring')
                    pl_string = '+'.join(questionable_hla_types)
                    hla_plstring_ques.text = pl_string
                    questionable_specificities_el.append(hla_plstring_ques)
                    
                    # Add ABL string
                    hla_ablstring_ques = create_namespaced_element('HLA-ablstring')
                    hla_ablstring_ques.text = create_abl_string(pl_string, conversion_map)
                    questionable_specificities_el.append(hla_ablstring_ques)
                interpretation_el.append(questionable_specificities_el)

                # Negative specificities
                negative_specificities_el = create_namespaced_element('negative-specificities')
                if negative_hla_types:
                    # Add PL string
                    hla_plstring_neg = create_namespaced_element('HLA-plstring')
                    pl_string = '+'.join(negative_hla_types)
                    hla_plstring_neg.text = pl_string
                    negative_specificities_el.append(hla_plstring_neg)
                    
                    # Add ABL string
                    hla_ablstring_neg = create_namespaced_element('HLA-ablstring')
                    hla_ablstring_neg.text = create_abl_string(pl_string, conversion_map)
                    negative_specificities_el.append(hla_ablstring_neg)
                interpretation_el.append(negative_specificities_el)
                
                # Remove existing interpretation element if present (schema allows 0 or 1)
                existing_interpretation = find_element_ns(assay, 'hamlns:interpretation')
                if existing_interpretation is not None:
                    assay.remove(existing_interpretation)
                
                # Add the new interpretation element to the assay
                assay.append(interpretation_el)
                
                # Apply proper indentation to the interpretation element
                indent_xml(interpretation_el)

    # Write the modified XML to the output file
    try:
        # Pretty print the entire tree
        indent_xml(root)
        tree.write(output_file_path, xml_declaration=True)
        print(f"Processed XML successfully written to {output_file_path}")
    except IOError as e:
        print(f"Error writing XML file: {e}")

if __name__ == "__main__":
    # Get input filename from command line argument or use default
    input_filename = sys.argv[1] if len(sys.argv) > 1 else 'HAML.xml'
    
    # Create output filename by inserting '_interpreted' before the extension
    base, ext = os.path.splitext(input_filename)
    output_filename = f"{base}_interpreted{ext}"
    
    print(f"Processing {input_filename}")
    process_xml_file(input_filename, output_filename)