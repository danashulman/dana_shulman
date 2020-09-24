import pandas as pd
import pickle
from sklearn import feature_selection
from urllib import request
import os


# from bs4 import BeautifulSoup
# from requests_html import HTMLSession

# Constants:
FIRST_WEB_PATH = "http://www.targetscan.org/cgi-bin/targetscan/vert_72/targetscan.cgi?species=Human&gid=&mir_sc" \
                 "=&mir_c=&mir_nc=&mir_vnc=&mirg={}"
TARGET_PREFIX = 'http://www.targetscan.org'
WEB_PREFIX = TARGET_PREFIX + '/cgi-bin/targetscan/vert_72/'
K = 200
URL_PREFIX = '<A HREF='
TARGETS_PATH = r'targets data'


if os.path.isfile("connection_problem.pkl"):
    connection_problem = pickle.load(open("connection_problem.pkl", 'rb'))
else:
    connection_problem = set()


def create_tag_dict():
    """
    Creates dictionary that holds the new ID of each patient and it's tag - sick or healthy
    :return: Saves the dictionary as a pickle file
    """
    df = pd.read_csv("miRNA/patients_info_sorted.csv")
    index_list = list(df['BUID'])
    patients_dict = {}
    for index in index_list:
        if int(df[df['BUID'] == index]['cogdx']) == 1:
            patients_dict[index] = 0
        else:
            patients_dict[index] = 1
    with open('patients_dict.pkl', 'wb') as f:
        pickle.dump(patients_dict, f)


def create_miRNA_dict():
    """
    Creates dictionary that holds the index of each miRNA
    :return: Saves the dictionary as a pickle file
    """
    miRNA_dict = {}
    miRNAs = list(pd.read_csv("miRNA/miRNA_combined_data.csv")['miRNA'])
    for i in range(len(miRNAs)):
        miRNA_dict[i] = miRNAs[i]
    with open('miRNA_dict.pkl', 'wb') as f:
        pickle.dump(miRNA_dict, f)


def extract_significant_features():
    """
    Creates a file that holds the most significant miRNAs
    """
    # Extract data and labels from the expression level of each miRNA table:
    df = pd.read_csv("miRNA/R code/man_hyp_normalized_counts.csv")
    old_patients_list = list(df.columns)[1:]
    patients_list = [int(s.strip('X')) for s in old_patients_list]
    df = df.iloc[:,1:].T
    patients_dict = pickle.load(open('patients_dict.pkl', 'rb'))
    labels = []
    for patient in patients_list:
        labels.append(patients_dict[patient])

    # Extracting most significant features (miRNAs):
    best_features = feature_selection.SelectKBest(score_func=feature_selection.f_classif, k='all')
    selector = best_features.fit(df, labels)
    df_scores = pd.DataFrame(selector.scores_)
    df_pvalues = pd.DataFrame(selector.pvalues_)
    df_columns = pd.DataFrame(df.columns)
    feature_scores = pd.concat([df_columns, df_scores, df_pvalues], axis=1)
    feature_scores.columns = ['miRNAs', 'Score', 'P-value']
    feature_scores = feature_scores.sort_values('Score')
    df = feature_scores.dropna()
    dict = pickle.load(open('miRNA_dict.pkl', 'rb'))
    df = df.replace({"miRNAs": dict})
    df = df.reindex(index=df.index[::-1])
    df = df[df['P-value'] <= 0.05]
    df.to_csv("most_significant_miRNAs.csv", index=False)


def download_miRNA_targets(miRNA_name):
    """
    Given a miRNA name, the table that holds the target genes of this miRNA is downloaded from TargetScan
    :param miRNA_name: the name of the wanted miRNA
    :return: downloads the table for the miRNA or prints an error message
    """
    url = FIRST_WEB_PATH.format(miRNA_name)
    try:
        response = request.urlopen(url)
        web_content = response.read()
        web_string = web_content.decode('ASCII')
        download_table_index = web_string.find('Download table', 0, len(web_string))
        url_begin_index = web_string.rindex(URL_PREFIX, 0, download_table_index - 1)
        url = WEB_PREFIX + web_string[url_begin_index + len(URL_PREFIX) + 1: download_table_index - 2]
        response = request.urlopen(url)
        web_content = response.read()
        web_string = web_content.decode('ASCII')
        xlsx_index = web_string.find('txt', 0, len(web_string))
        url_begin_index = web_string.find(URL_PREFIX, 0, xlsx_index)
        url = TARGET_PREFIX + web_string[url_begin_index + len(URL_PREFIX) + 1: xlsx_index + 3]
        request.urlretrieve(url, TARGETS_PATH + '\\{}'.format(miRNA_name + '.csv'))
    except Exception as ex:
        print(ex)
        print("couldn't reach the targets of the following miRNA: ", miRNA_name)
        connection_problem.add(miRNA_name)
        pickle.dump(connection_problem, open("connection_problem.pkl", 'wb'))


def corresponding_genes():
    """
    Combines and filters the targets and combines repetitions
    :return: creates a csv file for the corresponding genes
    """
    total = pd.DataFrame()
    for file_name in os.listdir(TARGETS_PATH):
        total = total.append(pd.read_csv(TARGETS_PATH + '\\' + file_name, sep='\t')[['Target gene',
                                                                                      'Representative transcript',
                                                                                      'Gene name',
                                                                                      'Cumulative weighted context++ score',
                                                                                      'Aggregate PCT']])

    total = total.dropna()
    total['Aggregate PCT'] = total['Aggregate PCT'].str.replace('<', '')
    total['Aggregate PCT'] = total['Aggregate PCT'].str.replace('>', '')
    total = total[total['Aggregate PCT'] != 'ORF']  # Think if I should keep them
    total['Aggregate PCT'] = total['Aggregate PCT'].astype(float)
    total['Cumulative weighted context++ score'] = total['Cumulative weighted context++ score'].astype(float)
    total = total[total['Cumulative weighted context++ score'] <= -0.2]  # decide the cutoff
    total = total[total['Aggregate PCT'] >= 0.8]  # decide the precision
    total = total[['Target gene', 'Representative transcript', 'Gene name']]
    # TODO: maybe filter by the name of the gene put of a list of associations
    total = total.groupby(['Target gene', 'Representative transcript', 'Gene name']).size().reset_index()
    total.columns = ['Target gene', 'Representative transcript', 'Gene name', 'Frequency']
    total = total.sort_values('Frequency', ascending=False)
    total.to_csv("corresponding_genes.csv", index=False)


# session = HTMLSession()
#
#
# def get_all_forms(url):
#     """Returns all form tags found on a web page's `url` """
#     # GET request
#     res = session.get(url)
#     # for javascript driven website
#     # res.html.render()
#     soup = BeautifulSoup(res.html.html, "html.parser")
#     return soup.find_all("form")
#
#
# def get_form_details(form):
#     """Returns the HTML details of a form,
#     including action, method and list of form controls (inputs, etc)"""
#     details = {}
#     # get the form action (requested URL)
#     action = form.attrs.get("action").lower()
#     # get the form method (POST, GET, DELETE, etc)
#     # if not specified, GET is the default in HTML
#     method = form.attrs.get("method", "get").lower()
#     # get all form inputs
#     inputs = []
#     for input_tag in form.find_all("input"):
#         # get type of input form control
#         input_type = input_tag.attrs.get("type", "text")
#         # get name attribute
#         input_name = input_tag.attrs.get("name")
#         # get the default value of that input tag
#         input_value =input_tag.attrs.get("value", "")
#         # add everything to that list
#         inputs.append({"type": input_type, "name": input_name, "value": input_value})
#     # put everything to the resulting dictionary
#     details["action"] = action
#     details["method"] = method
#     details["inputs"] = inputs
#     return details
#
#
# url = "https://string-db.org/"
# # get all form tags
# forms = get_all_forms(url)
# # iteratte over forms
# for i, form in enumerate(forms, start=1):
#     form_details = get_form_details(form)
#     print("="*50, f"form #{i}", "="*50)
#     print(form_details)

def create_data_for_deseq():
    """
    Creates the wanted table for the DESEQ2
    :return: sorted patients table that holds the RIN and samples' ID as well
    """
    df = pd.read_csv("../patients_info.csv").sort_values('BUID').iloc[:, 1:]
    df = df.reset_index(drop=True)
    old_cols = list(df.columns)
    new_cols = old_cols.copy()
    new_cols[0] = 'BUID'
    new_cols[1] = 'subjectID'
    df = df.reindex(columns=new_cols)
    rin_and_sample_ID = pd.read_csv("../more_patients_info.csv")[['RIN', 'sampleID']]
    combined = pd.concat([df, rin_and_sample_ID], axis=1)
    combined.to_csv("patients_info_sorted.csv", index=False)


def separate_miRNA_info_data():
    df_info = pd.read_csv("patients_info_sorted.csv")[['BUID', 'subjectID', 'miRNA_batch', 'brain_region',
                                                  'cogdx', 'age_death', 'msex', 'RIN', 'sampleID']]
    df_info = df_info.astype({'RIN': int, 'age_death': int})
    df_info = df_info[df_info['RIN'] >= 4]
    nuc_info = df_info[df_info['brain_region'] == 'NucAcc']
    hyp_info = df_info[df_info['brain_region'] == 'Hypothalamus']
    nuc_man_info = nuc_info[nuc_info['msex'] == 1]
    nuc_man_list = [str(x) for x in list(nuc_man_info['BUID'])]
    nuc_woman_info = nuc_info[nuc_info['msex'] == 0]
    nuc_woman_list = [str(x) for x in list(nuc_woman_info['BUID'])]
    hyp_man_info = hyp_info[hyp_info['msex'] == 1]
    hyp_man_list = [str(x) for x in list(hyp_man_info['BUID'])]
    hyp_woman_info = hyp_info[hyp_info['msex'] == 0]
    hyp_woman_list = [str(x) for x in list(hyp_woman_info['BUID'])]
    nuc_man_info.to_csv("nuc_man_info.csv", index=False)
    nuc_woman_info.to_csv("nuc_woman_info.csv", index=False)
    hyp_man_info.to_csv("hyp_man_info.csv", index=False)
    hyp_woman_info.to_csv("hyp_woman_info.csv", index=False)

    df_data = pd.read_csv("miRNA_combined_data.csv").iloc[:, 2:]
    nuc_man_data = df_data[nuc_man_list]
    nuc_woman_data = df_data[nuc_woman_list]
    hyp_man_data = df_data[hyp_man_list]
    hyp_woman_data = df_data[hyp_woman_list]
    nuc_man_data.to_csv("nuc_man_data.csv", index=False)
    nuc_woman_data.to_csv("nuc_woman_data.csv", index=False)
    hyp_man_data.to_csv("hyp_man_data.csv", index=False)
    hyp_woman_data.to_csv("hyp_woman_data.csv", index=False)


def change_names():
    main_dir_name = "EXAMPLE"
    for folder_name in os.listdir(main_dir_name):
        for file_name in os.listdir(main_dir_name + '/' + folder_name):
            if file_name.endswith(".csv"):
                os.rename(main_dir_name + '/' + folder_name + '/' + file_name, main_dir_name + '/' + folder_name + '/' + folder_name + '.csv')

change_names()
