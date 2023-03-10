import requests
import pp
import datetime

fullURL = ('https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cmass%2Ccc_sequence_caution%2Corganism_id%2Cannotation_score%2Ccc_miscellaneous%2Cprotein_existence%2Ccc_subunit%2Cgo_p%2Cgo%2Cgo_id%2Cstructure_3d%2Clit_pubmed_id%2Cdate_created%2Cdate_modified%2Cdate_sequence_modified%2Cxref_proteomes%2Cft_transmem%2Cft_topo_dom%2Cft_intramem%2Ccc_subcellular_location%2Cxref_ccds%2Cxref_embl%2Cxref_pir%2Cxref_refseq%2Cxref_alphafolddb%2Cxref_bmrb%2Cxref_pdb%2Cxref_pdbsum%2Cxref_sasbdb%2Cxref_pcddb%2Cxref_smr%2Cxref_biogrid%2Cxref_corum%2Cxref_complexportal%2Cxref_dip%2Cxref_elm%2Cxref_intact%2Cxref_mint%2Cxref_string%2Cxref_carbonyldb%2Cxref_depod%2Cxref_glyconnect%2Cxref_glygen%2Cxref_metosite%2Cxref_phosphositeplus%2Cxref_swisspalm%2Cxref_iptmnet%2Cxref_unicarbkb%2Cxref_interpro%2Cxref_genetree%2Cxref_cptac%2Cxref_epd%2Cxref_massive%2Cxref_maxqb%2Cxref_pride%2Cxref_paxdb%2Cxref_peptideatlas%2Cxref_promex%2Cxref_proteomicsdb%2Cxref_topdownproteomics%2Cxref_jpost&format=xlsx&query=%28taxonomy_id%3A2261%29')


result = requests.get(fullURL)

# save result separated by tab to xlsx with file containing current date, time, and organism name
with open('./Uniprot_input/uniprot_{}.xlsx'.format(datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")), 'wb') as f:
    f.write(result.content) # write content of result to file
    f.close() # close file

    