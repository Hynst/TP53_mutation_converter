import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.variantmapper
import hgvs.assemblymapper
import sys
import pandas as pd
import hgvs.validator
import hgvs.exceptions

##############################
#####     - TP53 -       #####
##### mutation converter #####
##############################

# Tool to convert variants anotations using HGVS standarts
# Two branches:
# 1) from coordinates to HGVS
#   run script as: python Convert_HGVS.py table1.txt "toHGVS"
#                  where table1.txt is tab delimited file with columns including:
#                  (Chr, Start, End, Ref, Alt, variant_table)
# 2) from HGVS_c to HGVS_g
#   run script as: python Convert_HGVS.py table1.txt "fromHGVS"
#                  where table1.txt is tab delimited file with columns including:
#                  (HGVS_c, variant_type)

# Load argument file
print("### Input file: " + str(sys.argv[1]))
print("### Convert to: " + str(sys.argv[2]))
print("### Run ID: " + str(sys.argv[3]))
variants_table = pd.read_table(sys.argv[1], sep="\t")

################################
### From coordinates to HGVS ###
################################
if str(sys.argv[2]) == "toHGVS":
    # NC_000017.10 NCBI ID for chromosome 17 in hg19
    assembly="NC_000017.10"
    START = list(variants_table['Start'])
    END = list(variants_table['End'])
    REF = list(variants_table['Ref'])
    ALT = list(variants_table['Alt'])

    # create hgvs from  Chrom, Pos, Ref, Alt
    final_table = pd.DataFrame(columns=['Chr','Start','End','Ref','Alt','HGVS_g','HGVS_c','HGVS_p'], index=range(len(START)))

    hp = hgvs.parser.Parser()
    hdp = hgvs.dataproviders.uta.connect()
    vr = hgvs.validator.Validator(hdp=hdp)
    # initialize the mapper for GRCh37 with splign-based alignments
    variantmapper = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37",alt_aln_method='splign',replace_reference=True)

    for var in xrange(len(START)):

        if variants_table['variant_type'][var] == "snv":

            # create HGVS_g
            hgvs_g = assembly + ":" + "g" + "." + str(START[var]) + REF[var] + ">" + ALT[var]
            print("Converting: " + str(hgvs_g))

            #hp = hgvs.parser.Parser()
            var_g = hp.parse_hgvs_variant(hgvs_g)

            transcripts = variantmapper.relevant_transcripts(var_g)
            var_c = variantmapper.g_to_c(var_g, 'NM_000546.5')
            var_p = variantmapper.c_to_p(var_c)

            # create final table
            final_table.loc[var].Chr = "17"
            final_table.loc[var].Start = START[var]
            final_table.loc[var].End = END[var]
            final_table.loc[var].Ref = REF[var]
            final_table.loc[var].Alt = ALT[var]
            final_table.loc[var].HGVS_g = var_g
            final_table.loc[var].HGVS_c = var_c
            final_table.loc[var].HGVS_p = var_p

        else:
            # create HGVS_g
            variant_type = str(variants_table['variant_type'][var])
            hgvs_g = assembly + ":" + "g" + "." + str(START[var]) + "_" + str(END[var]) + variant_type #+ str(REF[var])
            print("Converting: " + str(hgvs_g))

            #hp = hgvs.parser.Parser()
            var_g = hp.parse_hgvs_variant(hgvs_g)
            transcripts = variantmapper.relevant_transcripts(var_g)
            var_c = variantmapper.g_to_c(var_g, 'NM_000546.5')
            var_p = variantmapper.c_to_p(var_c)

            final_table.loc[var].Chr = "17"
            final_table.loc[var].Start = START[var]
            final_table.loc[var].End = END[var]
            final_table.loc[var].Ref = REF[var]
            final_table.loc[var].Alt = ALT[var]
            final_table.loc[var].HGVS_g = var_g
            final_table.loc[var].HGVS_c = var_c
            final_table.loc[var].HGVS_p = var_p

    final_table.to_csv(sys.argv[3] + "_toHGVS.txt", sep = "\t", index=False)
    print("### Final output file: " + str(sys.argv[3]) + "_toHGVS.txt")
########################
### HGVS_c to HGVS_g ###
########################
else:
    cDNA = list(variants_table['HGVS_c'])
    v_type = list(variants_table['variant_type'])
    final_table2 = pd.DataFrame(columns=['HGVS_c', 'HGVS_g', 'HGVS_p','Chr','Start','End', 'Ref', 'Alt'], index=range(len(cDNA)))

    hdp = hgvs.dataproviders.uta.connect()
    variantmapper = hgvs.variantmapper.VariantMapper(hdp)
    hp = hgvs.parser.Parser()

    # Map HGVS to genome
    for var2 in xrange(len(cDNA)):

        v_c = str(cDNA[var2])
        HGVS_c = hp.parse_hgvs_variant(v_c)
        print("Converting: " + str(HGVS_c))

        HGVS_g_str = variantmapper.c_to_g(HGVS_c, "NC_000017.10")
        HGVS_p_str = variantmapper.c_to_p(HGVS_c)

        ### convert SNV ###
        if str(v_type[var2]) == "snv":

            edit=str(HGVS_g_str.posedit.edit)

            # create final table
            final_table2.loc[var2].HGVS_c = v_c
            final_table2.loc[var2].HGVS_p = HGVS_p_str
            final_table2.loc[var2].HGVS_g = HGVS_g_str
            final_table2.loc[var2].Chr = "17"
            final_table2.loc[var2].Start = HGVS_g_str.posedit.pos
            final_table2.loc[var2].End = HGVS_g_str.posedit.pos
            final_table2.loc[var2].Ref = edit.split(">")[0]
            final_table2.loc[var2].Alt = edit.split(">")[1]

        ### convert INDEL ###
        else:
            pos=str(HGVS_g_str.posedit.pos)
            type=str(HGVS_g_str.posedit.edit)

            final_table2.loc[var2].HGVS_c = v_c
            final_table2.loc[var2].HGVS_g = HGVS_g_str
            final_table2.loc[var2].Chr = "17"
            if len(pos.split("_")) == 1:
                final_table2.loc[var2].Start = pos.split("_")[0]
                final_table2.loc[var2].End = pos.split("_")[0]
            else:
                final_table2.loc[var2].Start = pos.split("_")[1]
                final_table2.loc[var2].End = pos.split("_")[0]
            final_table2.loc[var2].Ref = "-"
            final_table2.loc[var2].Alt = "-"

    final_table2.to_csv(sys.argv[3] + "_fromHGVS.txt", sep = "\t", index=False)
    print("### Final output file: " + str(sys.argv[3]) + "_fromHGVS.txt")
