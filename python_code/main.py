import csv
from util import *


def set_ref_distr_to_CPOP(data_list, total_reduction, ets_share):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    total_reduction : float
        reduction target as decimal number (e.g. 55% = 0.55)
    ets_share : float
        share of emissions of EU ETS
    """
    df = data_list["Base_POP"][2018][1:] / data_list["Base_POP"][2018][1:].sum()
    data_list["Ref_distr"]["%-Share"][1:] = df
    data_list["Ref_distr"]["Mill t/y"] = (
        df
        * data_list["Base_C"][1990]["EU-27"]
        * (1 - total_reduction)
        * (1 - ets_share)
    )


def calc_goal_EU_2030(data_list, total_reduction, ets_share):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    total_reduction : float
        reduction target as decimal number (e.g. 55% = 0.55)
    ets_share : float
        share of emissions of EU ETS
    Returns
    -------
    float
        total carbon budget for EU  in 2030
    """
    budget = data_list["Base_C"][1990]["EU-27"] * (1 - total_reduction)
    return (budget - budget * ets_share) - LULUCF_EMISSIONS_1990_2030


def create_GDPPOP_distr_share_old(data_list, percentage, *args):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    percentage : float
        intensity factor to apply criterion
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    ----------
    """
    c_per_person_2018 = data_list["Base_GDP"][2019] / data_list["Base_POP"][2019]
    c_pop_ratio = c_per_person_2018 / c_per_person_2018["EU-27"]
    c_pop_distr = (
        data_list["Ref_distr"]["%-Share"] / c_pop_ratio.drop("EU-27")
    ).dropna()
    c_pop_distr = c_pop_distr / c_pop_distr.sum()
    c_pop_distr = (1 - percentage) * data_list["Ref_distr"][
        "%-Share"
    ] + percentage * c_pop_distr
    return c_pop_distr.dropna()


def create_alt_cap_indic_share(indicator: str, year: int):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    indicator : str
        name of governance indicator
    year : int
        year for indicator record
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    if indicator not in gov_indic_df.columns:
        print(f"indicator must be one of: {gov_indic_df.columns[3:]}")
        return

    df = gov_indic_eu_df[gov_indic_eu_df["year"] == year]
    eu_pop = data_list["Base_POP"][year]
    pop_share = eu_pop[1:] / eu_pop[0]
    pop_share.name = "pop_share"

    df = df.set_index("countryname")
    df = df.join(pop_share)
    indic_eu_avg = (df[indicator] * df["pop_share"]).sum()
    ratio = df[indicator] / indic_eu_avg
    ref_distr = data_list["Ref_distr"]["%-Share"]
    new_distr = ref_distr / ratio
    new_distr = new_distr.dropna() / new_distr.dropna().sum()
    return new_distr


def create_EU_GDPPOP_distr_share():
    """Calculates distribution for given criteria selection and formats to desired unit
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    eu_gdp = pd.read_excel("EU_GDP_approach_distribution.xlsx", index_col=1)
    eu_gdp = eu_gdp["EU_GDP_approach"].squeeze()
    sorted_countries = create_country_list(sorting_list)
    ref = (
        data_list["Base_C_ES"][2005].drop("EU-27")
        / data_list["Base_C_ES"][2005]["EU-27"]
    )
    eu_gdp_share = (ref * (1 - eu_gdp)).dropna()
    return (
        eu_gdp_share.sort_index(key=lambda x: x.map(sorted_countries)).dropna()
        / eu_gdp_share.sum()
    )


def create_min_welfare_distr_share(data_list, percentage, *args):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    percentage : float
        intensity factor to apply criterion
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    welf_distr = data_list["Base_min_welf"][1].drop("EU-27")
    data_list["Ref_distr"]["%-Share"].dropna(inplace=True)
    welf_distr = (1 - percentage) * data_list["Ref_distr"]["%-Share"] + (
        percentage * welf_distr
    )
    return welf_distr.dropna()


def create_renewables_distr_share(data_list, percentage, *args):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    percentage : float
        intensity factor to apply criterion
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    res_change = data_list["Base_RES"][2019].drop("EU-27") - data_list["Base_RES"][
        2005
    ].drop("EU-27")
    avg_wghtd_calc = res_change * data_list["Base_POP"][2019].drop("EU-27")
    avg_EU = avg_wghtd_calc.sum() / data_list["Base_POP"][2019].drop("EU-27").sum()
    diff_to_avg = res_change / avg_EU
    df = diff_to_avg * data_list["Ref_distr"]["%-Share"]
    return df.dropna() / df.sum()


def create_res_cap_distr_share(data_list, percentage, *args):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    percentage : float
        intensity factor to apply criterion
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    res_change = data_list["Base_RES"][2019].drop("EU-27") - data_list["Base_RES"][
        2005
    ].drop("EU-27")
    avg_wghtd_calc = res_change * data_list["Base_POP"][2019].drop("EU-27")
    avg_EU = avg_wghtd_calc.sum() / data_list["Base_POP"][2019].drop("EU-27").sum()
    diff_to_avg = res_change / avg_EU
    df = data_list["Ref_distr"]["%-Share"] / diff_to_avg
    return df.dropna() / df.sum()


def create_CPOP_distr_share(data_list, percentage, *args):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    percentage : float
        intensity factor to apply criterion
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    c_per_person_2018 = data_list["Base_C"][2019] / data_list["Base_POP"][2019]
    c_pop_ratio = c_per_person_2018 / c_per_person_2018["EU-27"]
    c_pop_distr = (
        data_list["Ref_distr"]["%-Share"] / c_pop_ratio.drop("EU-27")
    ).dropna()
    c_pop_distr = c_pop_distr / c_pop_distr.sum()
    c_pop_distr = (1 - percentage) * data_list["Ref_distr"][
        "%-Share"
    ] + percentage * c_pop_distr
    c_pop_distr.dropna(inplace=True)

    # Just calculate EPC shares
    df = data_list["Base_POP"][2019][1:] / data_list["Base_POP"][2019][1:].sum()
    return df


def create_CPOP_distr_share_Stefan(data_list, percentage, *args):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    percentage : float
        intensity factor to apply criterion
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    c_per_person_2019 = data_list["Base_C"][2019] / data_list["Base_POP"][2019]
    c_pop_ratio = c_per_person_2019 / c_per_person_2019["EU-27"]
    c_pop_distr = (
        data_list["Ref_distr"]["%-Share"] / c_pop_ratio.drop("EU-27")
    ).dropna()
    c_pop_distr = c_pop_distr / c_pop_distr.sum()
    c_pop_distr = (1 - percentage) * data_list["Ref_distr"][
        "%-Share"
    ] + percentage * c_pop_distr
    c_pop_distr.dropna(inplace=True)
    return c_pop_distr


def create_GDPPOP_distr_share(data_list, percentage, *args):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    percentage : float
        intensity factor to apply criterion
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    avg_wghtd_calc = data_list["Base_GDP"][2019].drop("EU-27")
    avg_EU = avg_wghtd_calc.sum() / data_list["Base_POP"][2019].drop("EU-27").sum()
    diff_to_avg = (
        data_list["Base_GDP"][2019].drop("EU-27")
        / data_list["Base_POP"][2019].drop("EU-27")
    ) / avg_EU
    df = (-diff_to_avg + 2) * data_list["Ref_distr"]["%-Share"]
    return df.dropna() / df.sum()


def create_CBUDGET_distr_share(data_list, percentage, goal_2030_EU):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    percentage : float
        intensity factor to apply criterion
    goal_2030_EU : float
        absolute carbon budget target for EU in 2030
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    years = 12
    c_budget = (data_list["Base_C_ES"][2019]["EU-27"] + goal_2030_EU) * years / 2
    c_pop_ratio = data_list["Base_POP"][2019] / data_list["Base_POP"][2019]["EU-27"]
    c_per_person_total = c_budget * c_pop_ratio
    c_per_person_total.drop("EU-27", inplace=True)
    c_2030 = c_per_person_total * 2 / years - data_list["Base_C_ES"][2019]
    c_2030 = c_2030 * goal_2030_EU / c_2030.sum()
    c_2030.dropna(inplace=True)

    neg_list = []
    for i in c_2030.items():
        if i[1] < 0:
            neg_list.append(i[0])

    remaining_emissions = -(c_2030[neg_list]).sum()
    c_2030 = c_2030 - (
        c_2030 * remaining_emissions / (c_2030.sum() + remaining_emissions)
    )
    c_2030[neg_list] *= 0
    c_pop_distr = c_2030 / goal_2030_EU
    c_pop_distr = (1 - percentage) * data_list["Ref_distr"][
        "%-Share"
    ] + percentage * c_pop_distr
    c_pop_distr.dropna(inplace=True)
    return c_pop_distr


def create_CBUDGET_1995_distr_share(
    data_list, percentage, goal_2030_EU, start_year=1995
):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    percentage : float
        intensity factor to apply criterion
    goal_2030_EU : float
        total carbon budget for EU in 2030
    start_year : int
        starting year to calculate already emitted carbon budget
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    years = 12
    c_budget = (data_list["Base_C_ES"][2019]["EU-27"] + goal_2030_EU) * years / 2
    c_budget -= data_list["Base_C_ES"][2019]["EU-27"]
    years_from_1990 = start_year - 1990
    df_cum_em = data_list["Base_C_ES_1990"][
        data_list["Base_C_ES_1990"].columns[years_from_1990:]
    ].sum(axis=1)
    c_budget = c_budget + df_cum_em.drop("EU-27").sum()
    c_pop_ratio = data_list["Base_POP"][2019] / data_list["Base_POP"][2019]["EU-27"]
    c_per_person_total = c_budget * c_pop_ratio - df_cum_em.drop("EU-27")
    c_per_person_total.drop("EU-27", inplace=True)
    c_per_person_total.dropna(inplace=True)

    c_2030 = (
        c_per_person_total + data_list["Base_C_ES"][2019].dropna().drop("EU-27")
    ) * 2 / years - data_list["Base_C_ES"][2019]
    c_2030 = c_2030 * goal_2030_EU / c_2030.sum()
    c_2030.dropna(inplace=True)

    neg_list = []
    for i in c_2030.items():
        if i[1] < 0:
            neg_list.append(i[0])

    remaining_emissions = -(c_2030[neg_list]).sum()
    c_2030 = c_2030 - (
        c_2030 * remaining_emissions / (c_2030.sum() + remaining_emissions)
    )
    c_2030[neg_list] *= 0
    c_pop_distr = c_2030 / goal_2030_EU
    c_pop_distr = (1 - percentage) * data_list["Ref_distr"][
        "%-Share"
    ] + percentage * c_pop_distr
    c_pop_distr.dropna(inplace=True)
    return c_pop_distr


def create_CBUDGET_1995_Stefan_distr_share(data_list):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    data_list["Base_C"].T.sum()
    budget_per_pop = (
        data_list["Base_C"].T.loc[1995:, :].sum()
        / data_list["Base_POP"].T.loc[1995:2019, :].sum()
    )
    diff_vector = budget_per_pop / budget_per_pop["EU-27"]
    shares_2030 = data_list["Ref_distr"]["%-Share"] / diff_vector
    shares_2030.drop("EU-27", inplace=True)
    return shares_2030 / shares_2030.sum()


def create_CBUDGET_2005_distr_share(data_list, percentage, goal_2030_EU):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    percentage : float
        intensity factor to apply criterion
    goal_2030_EU : float
        absolute carbon budget target for EU in 2030
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    years = 12
    c_budget = (data_list["Base_C_ES"][2019]["EU-27"] + goal_2030_EU) * years / 2
    c_budget -= data_list["Base_C_ES"][2019]["EU-27"]
    c_budget = (
        c_budget
        + (
            data_list["Base_1995"][2019].drop("EU-27")
            - data_list["Base_1995"][2004].drop("EU-27")
        ).sum()
    )
    c_pop_ratio = data_list["Base_POP"][2019] / data_list["Base_POP"][2019]["EU-27"]
    c_per_person_total = c_budget * c_pop_ratio - (
        data_list["Base_1995"][2019].drop("EU-27")
        - data_list["Base_1995"][2004].drop("EU-27")
    )
    c_per_person_total.drop("EU-27", inplace=True)
    c_per_person_total.dropna(inplace=True)

    c_2030 = (
        c_per_person_total + data_list["Base_C_ES"][2019].dropna().drop("EU-27")
    ) * 2 / years - data_list["Base_C_ES"][2019]
    c_2030 = c_2030 * goal_2030_EU / c_2030.sum()
    c_2030.dropna(inplace=True)

    neg_list = []
    for i in c_2030.items():
        if i[1] < 0:
            neg_list.append(i[0])

    remaining_emissions = -(c_2030[neg_list]).sum()
    c_2030 = c_2030 - (
        c_2030 * remaining_emissions / (c_2030.sum() + remaining_emissions)
    )
    c_2030[neg_list] *= 0
    c_pop_distr = c_2030 / goal_2030_EU
    c_pop_distr = (1 - percentage) * data_list["Ref_distr"][
        "%-Share"
    ] + percentage * c_pop_distr
    c_pop_distr.dropna(inplace=True)
    return c_pop_distr


def create_inh_benefits_distr_share(data_list, percentage, goal_2030_EU):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    percentage : float
        intensity factor to apply criterion
    goal_2030_EU : float
        absolute carbon budget target for EU in 2030
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    years = 12
    c_budget = (data_list["Base_C_ES"][2019]["EU-27"] + goal_2030_EU) * years / 2
    c_budget -= data_list["Base_C_ES"][2019]["EU-27"]
    c_budget = c_budget + data_list["Base_benefits"]["Mill. Tonnes"]["EU-27"]
    c_pop_ratio = data_list["Base_POP"][2019] / data_list["Base_POP"][2019]["EU-27"]
    c_per_person_total = c_budget * c_pop_ratio - (
        data_list["Base_benefits"]["Mill. Tonnes"]
    )
    c_per_person_total.drop("EU-27", inplace=True)
    c_per_person_total.dropna(inplace=True)
    c_2030 = (
        c_per_person_total + data_list["Base_C_ES"][2019].dropna().drop("EU-27")
    ) * 2 / years - data_list["Base_C_ES"][2019]
    c_2030 = c_2030 * goal_2030_EU / c_2030.sum()
    c_2030.dropna(inplace=True)

    neg_list = []
    for i in c_2030.items():
        if i[1] < 0:
            neg_list.append(i[0])

    remaining_emissions = -(c_2030[neg_list]).sum()
    c_2030 = c_2030 - (
        c_2030 * remaining_emissions / (c_2030.sum() + remaining_emissions)
    )
    c_2030[neg_list] *= 0
    c_pop_distr = c_2030 / goal_2030_EU
    c_pop_distr = (1 - percentage) * data_list["Ref_distr"][
        "%-Share"
    ] + percentage * c_pop_distr
    c_pop_distr.dropna(inplace=True)
    return c_pop_distr


def zero_restriction_correction_loop(data_list, goal_EU_2030, qual_distr):
    """Applies the zero restriction mechanism, to ensure that the selected distribution
    does not allow for positive emission changes for 2030 compared to 2005. Remaining
    emissions are redistributed.

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    goal_2030_EU : float
        absolute carbon budget target for EU in 2030
    qual_distr : pd.DataFrame
        DataFrame with the budget distribution over EU MS
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    loop_bool = True
    change_2005_2030 = (qual_distr * goal_EU_2030) / data_list["Base_C_ES"][2005] - 1
    while loop_bool:
        zero_list = []
        for i in change_2005_2030.items():
            if i[1] >= 0:
                zero_list.append(i[0])

        qual_2030_absolute = qual_distr * goal_EU_2030

        remaining_emissions = (
            qual_2030_absolute[zero_list] - data_list["Base_C_ES"][2005][zero_list]
        ).sum()
        qual_2030_absolute = qual_2030_absolute + (
            qual_2030_absolute
            * remaining_emissions
            / (
                qual_2030_absolute.sum()
                - data_list["Base_C_ES"][2005][zero_list].sum()
                - remaining_emissions
            )
        )
        qual_2030_absolute[zero_list] = data_list["Base_C_ES"][2005][zero_list]
        change_2005_2030 = qual_2030_absolute / data_list["Base_C_ES"][2005] - 1
        if change_2005_2030.max() <= 0:
            loop_bool = False

    qual_distr = qual_2030_absolute / goal_EU_2030
    return qual_distr


def calc_qualified_combo_new(qual_par_list, data_list, goal_EU_2030, **kwargs):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    data_list : pd.DataFrame
        emission country data
    goal_EU_2030 : float
        equal per capita budget criterion
    Returns
    -------
    pd.DataFrame
        DataFrame with relative emission budget distribution
    """
    kwargs.setdefault("zero_res", True)

    if all_zero(qual_par_list.values()):
        cpop = create_CPOP_distr_share(data_list, qual_par_list["CPOP"])
        cpop = data_list["Ref_distr"]["%-Share"]
        if "zero_res" in kwargs.keys():
            if kwargs.get("zero_res") == True:
                cpop = zero_restriction_correction_loop(data_list, goal_EU_2030, cpop)
        return cpop

    qual_distr = []

    if "year_sel" in kwargs.keys():
        year = kwargs.get("year_sel")
    else:
        year = 2020

    if kwargs.get("zero_res") == False:

        if qual_par_list["CPOP"] != 0:
            cpop = create_CPOP_distr_share(data_list, 1)
            qual_distr.append(cpop * qual_par_list["CPOP"])

        if qual_par_list["EU_GDPPOP"] != 0:
            eu_gdppop = create_EU_GDPPOP_distr_share()
            qual_distr.append(eu_gdppop * qual_par_list["EU_GDPPOP"])

        if qual_par_list["CPOP_stefan"] != 0:
            eu_gdppop = create_CPOP_distr_share_Stefan(data_list, 1)
            qual_distr.append(eu_gdppop * qual_par_list["CPOP_stefan"])

        if qual_par_list["CBUD_stefan"] != 0:
            eu_gdppop = create_CBUDGET_1995_Stefan_distr_share(data_list)
            qual_distr.append(eu_gdppop * qual_par_list["CBUD_stefan"])

        if qual_par_list["GDPPOP_div"] != 0:
            gdppop = create_GDPPOP_distr_share_old(data_list, 1)
            qual_distr.append(gdppop * qual_par_list["GDPPOP_div"])

        if qual_par_list["GDPPOP"] != 0:
            gdppop = create_GDPPOP_distr_share(data_list, 1)
            qual_distr.append(gdppop * qual_par_list["GDPPOP"])

        if qual_par_list["MIN_WELF"] != 0:
            welf = create_min_welfare_distr_share(data_list, 1)
            qual_distr.append(welf * qual_par_list["MIN_WELF"])

        if qual_par_list["RES"] != 0:
            res = create_renewables_distr_share(data_list, 1)
            qual_distr.append(res * qual_par_list["RES"])

        if qual_par_list["CBUDGET"] != 0:
            cbudget = create_CBUDGET_distr_share(data_list, 1, goal_EU_2030)
            qual_distr.append(cbudget * qual_par_list["CBUDGET"])

        if qual_par_list["1995"] != 0:
            cbudget_1995 = create_CBUDGET_1995_distr_share(data_list, 1, goal_EU_2030)
            qual_distr.append(cbudget_1995 * qual_par_list["1995"])

        if qual_par_list["2005"] != 0:
            cbudget_2005 = create_CBUDGET_2005_distr_share(data_list, 1, goal_EU_2030)
            qual_distr.append(cbudget_2005 * qual_par_list["2005"])

        if qual_par_list["BENEFITS"] != 0:
            benefits = create_inh_benefits_distr_share(data_list, 1, goal_EU_2030)
            qual_distr.append(benefits * qual_par_list["BENEFITS"])

        if qual_par_list["RES_cap"] != 0:
            res = create_res_cap_distr_share(data_list, 1)
            qual_distr.append(res * qual_par_list["RES_cap"])

        alt_indic_list = [
            i for i in qual_par_list.keys() if i in gov_indic_eu_df.columns
        ]
        if alt_indic_list:
            for i in alt_indic_list:
                result = create_alt_cap_indic_share(i, year)
                qual_distr.append(result * qual_par_list[i])

    else:
        if qual_par_list["CPOP"] != 0:
            cpop = create_CPOP_distr_share(data_list, 1)
            cpop = zero_restriction_correction_loop(data_list, goal_EU_2030, cpop)
            qual_distr.append(cpop * qual_par_list["CPOP"])

        if qual_par_list["CPOP_stefan"] != 0:
            eu_gdppop = zero_restriction_correction_loop(
                data_list, goal_EU_2030, create_CPOP_distr_share_Stefan(data_list, 1)
            )
            qual_distr.append(eu_gdppop * qual_par_list["CPOP_stefan"])

        if qual_par_list["CBUD_stefan"] != 0:
            eu_gdppop = zero_restriction_correction_loop(
                data_list,
                goal_EU_2030,
                create_CBUDGET_1995_Stefan_distr_share(data_list),
            )
            qual_distr.append(eu_gdppop * qual_par_list["CBUD_stefan"])

        if qual_par_list["GDPPOP"] != 0:
            gdppop = zero_restriction_correction_loop(
                data_list, goal_EU_2030, create_GDPPOP_distr_share(data_list, 1)
            )
            qual_distr.append(gdppop * qual_par_list["GDPPOP"])

        if qual_par_list["EU_GDPPOP"] != 0:
            eu_gdppop = create_EU_GDPPOP_distr_share()
            qual_distr.append(eu_gdppop * qual_par_list["EU_GDPPOP"])

        if qual_par_list["GDPPOP_div"] != 0:
            gdppop = zero_restriction_correction_loop(
                data_list, goal_EU_2030, create_GDPPOP_distr_share_old(data_list, 1)
            )
            qual_distr.append(gdppop * qual_par_list["GDPPOP_div"])

        if qual_par_list["MIN_WELF"] != 0:
            welf = zero_restriction_correction_loop(
                data_list, goal_EU_2030, create_min_welfare_distr_share(data_list, 1)
            )
            qual_distr.append(welf * qual_par_list["MIN_WELF"])

        if qual_par_list["RES"] != 0:
            res = zero_restriction_correction_loop(
                data_list, goal_EU_2030, create_renewables_distr_share(data_list, 1)
            )
            qual_distr.append(res * qual_par_list["RES"])

        if qual_par_list["CBUDGET"] != 0:
            cbudget = zero_restriction_correction_loop(
                data_list,
                goal_EU_2030,
                create_CBUDGET_distr_share(data_list, 1, goal_EU_2030),
            )
            qual_distr.append(cbudget * qual_par_list["CBUDGET"])

        if qual_par_list["1995"] != 0:
            cbudget_1995 = zero_restriction_correction_loop(
                data_list,
                goal_EU_2030,
                create_CBUDGET_1995_distr_share(data_list, 1, goal_EU_2030),
            )
            qual_distr.append(cbudget_1995 * qual_par_list["1995"])

        if qual_par_list["2005"] != 0:
            cbudget_2005 = zero_restriction_correction_loop(
                data_list,
                goal_EU_2030,
                create_CBUDGET_2005_distr_share(data_list, 1, goal_EU_2030),
            )
            qual_distr.append(cbudget_2005 * qual_par_list["2005"])

        if qual_par_list["BENEFITS"] != 0:
            benefits = zero_restriction_correction_loop(
                data_list,
                goal_EU_2030,
                create_inh_benefits_distr_share(data_list, 1, goal_EU_2030),
            )
            qual_distr.append(benefits * qual_par_list["BENEFITS"])

        if qual_par_list["RES_cap"] != 0:
            res = zero_restriction_correction_loop(
                data_list, goal_EU_2030, create_res_cap_distr_share(data_list, 1)
            )
            qual_distr.append(res * qual_par_list["RES_cap"])

        alt_indic_list = [
            i for i in qual_par_list.keys() if i in gov_indic_eu_df.columns
        ]
        if alt_indic_list:
            for i in alt_indic_list:
                result = zero_restriction_correction_loop(
                    data_list, goal_EU_2030, create_alt_cap_indic_share(i, year)
                )
                qual_distr.append(result * qual_par_list[i])

    total_distr = qual_distr[0].copy()
    for i in qual_distr:
        total_distr += i
    return total_distr - qual_distr[0]


def calculate_distribution(
    cpop,
    res,
    cbud,
    gdp,
    welf,
    cbud_95,
    cbud_05,
    inh_benef,
    eu_gdp,
    gdp_div,
    cbud_stefan,
    cpop_stefan,
    gee_n,
    res_cap,
    **kwargs
):
    """Calculates distribution for given criteria selection and formats to desired unit

    Parameters
    ----------
    cpop : float
        equal per capita budget criterion
    res : float
        equal per capita budget criterion
    cbud : float
        equal per capita budget criterion
    gdp : float
        gdp per capita criterion
    welf: float
        basic needs criterion
    cbud_95 : float
        historical emissions 1995 criterion
    cbud_05 : float
        historical emissions 2005 criterion
    inh_benef : float
        inherited benefits of emissions criterion
    eu_gdp : float
        EU implementation approximation
    gdp_div : float
        GDP per capita criterion
    cbud_stefan : float
        equal carbon budget criterion
    cpop_stefan : float
        equal per capita criterion "old"
    gee_n : float
        government effectiveness criterion
    res_cap : float
        renewable growth capacity criterion
    key provision : bool
        apply provision mechanism
    **kwargs
        "output" unit can be set to "reduction", "share" or "absolute",
        "provision"=True to activate provision mechanism,
        "goal_EU_2030" absolute emission target for EU 2030
    """
    kwargs.setdefault("provision", True)
    output = kwargs.get("output")
    goal_EU_2030 = kwargs.get("goal_EU_2030")
    qual_dict = {
        "CPOP": cpop,
        "CPOP_stefan": cpop_stefan,
        "RES": res,
        "CBUDGET": cbud,
        "GDPPOP": gdp,
        "EU_GDPPOP": eu_gdp,
        "GDPPOP_div": gdp_div,
        "MIN_WELF": welf,
        "1995": cbud_95,
        "CBUD_stefan": cbud_stefan,
        "2005": cbud_05,
        "BENEFITS": inh_benef,
        "gee_n": gee_n,
        "RES_cap": res_cap,
    }
    distr_combined = calc_qualified_combo_new(qual_dict, data_list, **kwargs)
    if "EU-27" in distr_combined:
        distr_combined.drop("EU-27", inplace=True)
    distr_combined *= goal_EU_2030
    if output == "reduction":
        distr_combined = distr_combined / data_list["Base_C_ES"][2005].dropna().drop("EU-27") - 1
    if output == "share":
        distr_combined = distr_combined / goal_EU_2030
    distr_combined.index.name = None

    sorted_countries = create_country_list(sorting_list)
    return distr_combined.sort_index(key=lambda x: x.map(sorted_countries))


def create_results_csv(**kwargs):
    """Function to create data for ternary graphs with 3 selectable interpretations
    """
    start = kwargs.get("start")
    stop = kwargs.get("stop")
    step = kwargs.get("step")
    qual_dict = kwargs.get("qual_dict")
    qual_var_index = map_qual_to_index(kwargs.get("qual_var"))
    qual_var_index1 = map_qual_to_index(kwargs.get("qual_var1"))
    qual_var_index2 = map_qual_to_index(kwargs.get("qual_var2"))
    file_name = kwargs.get("file_name") + ".csv"

    array = [
        qual_dict["CPOP"],
        qual_dict["RES"],
        qual_dict["CBUDGET"],
        qual_dict["GDPPOP"],
        qual_dict["MIN_WELF"],
        qual_dict["1995"],
        qual_dict["2005"],
        qual_dict["BENEFITS"],
        qual_dict["EU_GDPPOP"],
        qual_dict["GDPPOP_div"],
        qual_dict["CBUD_stefan"],
        qual_dict["CPOP_stefan"],
        qual_dict["gee_n"],
        qual_dict["RES_cap"],
    ]

    str_buf = ""
    row_list = []
    header = [
        "Target reduction",
        "ETS share",
        "CPOP",
        "GDPPOP",
        "RES",
        "CBUDGET",
        "MIN_WELF",
        "1995",
        "2005",
        "BENEFITS",
        "EU_GDPPOP",
        "GDPPOP_div",
        "CBUD_stefan",
        "CPOP_stefan",
        "gee_n",
        "RES_cap",
    ]
    header = ["countries", "value"] + header

    with open(file_name, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(header)

    par_val = start
    par_val1 = start
    par_val2 = start

    while par_val2 < stop:
        array[qual_var_index2] = par_val2
        while par_val1 < stop:
            array[qual_var_index1] = par_val1
            str_buf = ""
            while par_val < stop:
                row_list = []
                array[qual_var_index] = par_val

                x = calculate_distribution(*array, **kwargs)

                df_temp = pd.DataFrame(x)
                df_temp.reset_index(inplace=True)
                df_temp.rename(columns={"index": "countries", 0: "value"}, inplace=True)
                df_temp.countries[4] = "Germany"
                df_temp.countries = df_temp.countries.apply(
                    lambda x: pycountry.countries.lookup(x).alpha_3
                )
                df_temp.set_index("countries", inplace=True)

                df_temp["Target reduction"] = kwargs.get("reduction")
                df_temp["ETS share"] = kwargs.get("ets_share")
                df_temp["CPOP"] = array[0]
                df_temp["GDPPOP"] = array[3]
                df_temp["RES"] = array[1]
                df_temp["CBUDGET"] = array[2]
                df_temp["MIN_WELF"] = array[4]
                df_temp["1995"] = array[5]
                df_temp["2005"] = array[6]
                df_temp["BENEFITS"] = array[7]
                df_temp["EU_GDPPOP"] = array[8]
                df_temp["GDPPOP_div"] = array[9]
                df_temp["CBUD_stefan"] = array[10]
                df_temp["CPOP_stefan"] = array[11]
                df_temp["gee_n"] = array[12]
                df_temp["RES_cap"] = array[13]

                str_buf = df_temp.to_csv(header=False)

                for i in str_buf.split("\r\n"):
                    row_list.append([i.split(",")])
                row_list = row_list[:-1]
                with open(file_name, "a", newline="") as file:
                    writer = csv.writer(file)
                    writer.writerows(i[0] for i in row_list)

                par_val += step
            par_val = start
            par_val1 += step
        par_val1 = start
        par_val2 += step
