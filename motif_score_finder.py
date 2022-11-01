import pandas as pd
sites_table = pd.read_csv("PTM_Table.csv", low_memory=False)
dict_sites_PTM = {}
for row, column in sites_table.iterrows():
    UID = column[2]
    site_pos = int(column[0].split("-")[1])
    if UID in dict_sites_PTM:
        dict_sites_PTM[UID].append(site_pos)
    if UID not in dict_sites_PTM:
        dict_sites_PTM[UID] = [site_pos]



interpro_result = pd.read_table("interproresults_withgo_pa.tsv", sep = "\t")
#interpro_result.columns = ["Protein_accession", "Sequence_MD5_Digest", "Length", "Analysis", "Signature_accession", "Description", "Start_pos", "End_pos", "E-Val", "Status", "Date", "Interpro_anno", "Interpro_anno_2"]
interpro_result_pos = {}
for row, column in interpro_result.iterrows():
    UID = column[0].split("|")[1]
    start = int(column[6])
    stop = int(column[7])
    eval = column[8]
    pos_list = (start, stop, eval)
    if UID in interpro_result_pos:
        interpro_result_pos[UID].append(pos_list)
    if UID not in interpro_result_pos:
        interpro_result_pos[UID] = [pos_list]


result_appended = []
for UID_PTM in dict_sites_PTM:
    if UID in interpro_result_pos:
        list_pos = dict_sites_PTM[UID]
        for pos in list_pos:
            for result_pos in interpro_result_pos[UID]:
                if pos in range(result_pos[0], result_pos[1]):
                    result = []
                    result.append(UID)
                    result.append(pos)
                    result.append(result_pos[0])
                    result.append(result_pos[1])
                    result.append(result_pos[2])
                    result_appended.append(result)
result_df = pd.DataFrame(result_appended)
result_df.columns = ["UID","Site_position","Motif_start_pos", "Motif_stop_pos","E-value"]
result_df.to_csv("Results_with_GO_PA")