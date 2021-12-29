from __future__ import absolute_import
from __future__ import print_function

import argparse
import yaml

from mimic3benchmark.mimic3csv import *
from mimic3benchmark.preprocessing import add_hcup_ccs_2015_groups, make_phenotype_label_matrix
from mimic3benchmark.util import dataframe_from_csv

parser = argparse.ArgumentParser(description='Extract per-subject data from MIMIC-III CSV files.')
parser.add_argument('mimic3_path', type=str, help='Directory containing MIMIC-III CSV files.')
# 含有mimic-iii csv文件的路径
parser.add_argument('output_path', type=str, help='Directory where per-subject data should be written.')
# 输出路径
parser.add_argument('--event_tables', '-e', type=str, nargs='+', help='Tables from which to read events.',
                    default=['CHARTEVENTS', 'LABEVENTS', 'OUTPUTEVENTS'])
parser.add_argument('--phenotype_definitions', '-p', type=str,
                    default=os.path.join(os.path.dirname(__file__), '../resources/hcup_ccs_2015_definitions.yaml'),
                    help='YAML file with phenotype definitions.')
parser.add_argument('--itemids_file', '-i', type=str, help='CSV containing list of ITEMIDs to keep.')
parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', help='Verbosity in output')
parser.add_argument('--quiet', '-q', dest='verbose', action='store_false', help='Suspend printing of details')
parser.set_defaults(verbose=True)
parser.add_argument('--test', action='store_true', help='TEST MODE: process only 1000 subjects, 1000000 events.')
args, _ = parser.parse_known_args()

if not os.path.exists(args.output_path):
    os.makedirs(args.output_path)

patients = read_patients_table(args.mimic3_path)
# 病人记录
admits = read_admissions_table(args.mimic3_path)
# 医院访问记录
stays = read_icustays_table(args.mimic3_path)
# ICU访问记录
if args.verbose:
# 是否设置为详细输出
    print('START:\n\tICUSTAY_IDs: {}\n\tHADM_IDs: {}\n\tSUBJECT_IDs: {}'.format(stays.ICUSTAY_ID.unique().shape[0],
          stays.HADM_ID.unique().shape[0], stays.SUBJECT_ID.unique().shape[0]))

stays = remove_icustays_with_transfers(stays)
# 移除有过病房转移的ICU记录
if args.verbose:
    print('REMOVE ICU TRANSFERS:\n\tICUSTAY_IDs: {}\n\tHADM_IDs: {}\n\tSUBJECT_IDs: {}'.format(stays.ICUSTAY_ID.unique().shape[0],
          stays.HADM_ID.unique().shape[0], stays.SUBJECT_ID.unique().shape[0]))

stays = merge_on_subject_admission(stays, admits)
stays = merge_on_subject(stays, patients)
# 病人、医院访问记录和ICU记录 INNER JOIN
stays = filter_admissions_on_nb_icustays(stays)
# 筛选出ICU访问记录在指定范围内的医院访问记录
if args.verbose:
    print('REMOVE MULTIPLE STAYS PER ADMIT:\n\tICUSTAY_IDs: {}\n\tHADM_IDs: {}\n\tSUBJECT_IDs: {}'.format(stays.ICUSTAY_ID.unique().shape[0],
          stays.HADM_ID.unique().shape[0], stays.SUBJECT_ID.unique().shape[0]))

stays = add_age_to_icustays(stays)
# 计算进入ICU时的年龄
stays = add_inunit_mortality_to_icustays(stays)
# 添加是否在ICU内死亡这一列
stays = add_inhospital_mortality_to_icustays(stays)
# 添加是否在医院内死亡这一列
stays = filter_icustays_on_age(stays)
# 年龄筛选（已进入ICU时的年龄为准）
if args.verbose:
    print('REMOVE PATIENTS AGE < 18:\n\tICUSTAY_IDs: {}\n\tHADM_IDs: {}\n\tSUBJECT_IDs: {}'.format(stays.ICUSTAY_ID.unique().shape[0],
          stays.HADM_ID.unique().shape[0], stays.SUBJECT_ID.unique().shape[0]))

stays.to_csv(os.path.join(args.output_path, 'all_stays.csv'), index=False)
diagnoses = read_icd_diagnoses_table(args.mimic3_path)
# 入院记录对应的诊断信息
diagnoses = filter_diagnoses_on_stays(diagnoses, stays)
# 每一次ICU记录对应的诊断信息
diagnoses.to_csv(os.path.join(args.output_path, 'all_diagnoses.csv'), index=False)
count_icd_codes(diagnoses, output_path=os.path.join(args.output_path, 'diagnosis_counts.csv'))
# 每一种病对应的ICU记录数

phenotypes = add_hcup_ccs_2015_groups(diagnoses, yaml.load(open(args.phenotype_definitions, 'r'), yaml.FullLoader))
make_phenotype_label_matrix(phenotypes, stays).to_csv(os.path.join(args.output_path, 'phenotype_labels.csv'),
                                                      index=False, quoting=csv.QUOTE_NONNUMERIC)

if args.test:
    pat_idx = np.random.choice(patients.shape[0], size=1000)
    patients = patients.iloc[pat_idx]
    stays = stays.merge(patients[['SUBJECT_ID']], left_on='SUBJECT_ID', right_on='SUBJECT_ID')
    args.event_tables = [args.event_tables[0]]
    print('Using only', stays.shape[0], 'stays and only', args.event_tables[0], 'table')

subjects = stays.SUBJECT_ID.unique()
break_up_stays_by_subject(stays, args.output_path, subjects=subjects)
break_up_diagnoses_by_subject(phenotypes, args.output_path, subjects=subjects)
items_to_keep = set(
    [int(itemid) for itemid in dataframe_from_csv(args.itemids_file)['ITEMID'].unique()]) if args.itemids_file else None
for table in args.event_tables:
    read_events_table_and_break_up_by_subject(args.mimic3_path, table, args.output_path, items_to_keep=items_to_keep,
                                              subjects_to_keep=subjects)
