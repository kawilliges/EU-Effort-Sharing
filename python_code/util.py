import scipy.stats as scp
import pandas as pd
import pycountry
from countrygroups import EUROPEAN_UNION


data_list = pd.read_excel(
    "EU_data_new.xlsx", skiprows=[0, 2, 4], sheet_name=None, engine="openpyxl"
)
data_list["Base_C"].rename(columns={"GHG Emissions Total": "countries"}, inplace=True)
sorting_list = list(data_list["Base_C"]["countries"].items())

data_list = pd.read_excel(
    "EU_data_new.xlsx",
    skiprows=[0, 2, 4],
    sheet_name=None,
    index_col=0,
    engine="openpyxl",
)

gov_indic_df = pd.read_csv("capacity_indicators_18_01_2022.csv")
gov_indic_eu_df = gov_indic_df[
    (gov_indic_df["iso3c"] == "ROM") | (gov_indic_df["iso3c"].isin(EUROPEAN_UNION))
]
gov_indic_eu_df = gov_indic_eu_df[gov_indic_eu_df["icrg_qog"].notna()]
gov_indic_eu_df = gov_indic_eu_df.replace(
    {"Germany ": "Germany", "Slovak Republic": "Slovakia", "Czech Republic": "Czechia"}
)

ref_goal_EU_2030 = data_list["Ref_distr"]["Mill t/y"]
ref_goal_EU_2030_list = ref_goal_EU_2030[1:]
ref_goal_EU_2030 = ref_goal_EU_2030.sum()

LULUCF_EMISSIONS_1990_2030 = -98.76467000


def all_zero(iterable):
    """Helper to determine if all values in iterable are 0

    Parameters
    ----------
    iterable : list
        iterable object
    Returns
    -------
    bool
        True if all zero, else False
    """
    iterable = iter(iterable)
    try:
        first = 0
    except StopIteration:
        return True
    return all(first == x for x in iterable)


def create_country_list(sorting_list):
    """Creates dictionary of country order in a list

    Parameters
    ----------
    sorting_list : list
        list with countries
    Returns
    -------
    dict
        dictionary with position as key and country name as value
    """
    key = []
    value = []
    for i, j in sorting_list:
        key.append(j)
        value.append(i)
    return (dict(zip(key, value)))


def map_qual_to_index(qual):
    """Helper to get index of interpretation in the standardized interpretation array

    Parameters
    ----------
    qual : str
        interpretation name (e.g. "CPOP")
    Returns
    -------
    int
        index of interpretation
    """
    map_dict = {
        "CPOP": 1,
        "RES": 2,
        "CBUDGET": 3,
        "GDPPOP": 4,
        "MIN_WELF": 5,
        "1995": 6,
        "2005": 7,
        "BENEFITS": 8,
        "EU_GDPPOP": 9,
        "GDPPOP_div": 10,
        "CBUD_stefan" : 11,
        "CPOP_stefan" : 12,
        "gee_n" : 13,
        "RES_cap": 14
    }
    return map_dict[qual] - 1

