#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" trifid/visualization/figures.py
https://gitlab.com/fpozoc/trifid/-/blob/master/trifid/models/utils.py

TRIFID useful figures.

Classes and functions:
    * cat_appris_order
    * cat_transcript_type
    * config_altair
    * create_categories
    * explain_prediction
    * plot_appris_histogram
    * plot_feature_importances
    * plot_learning_curve
    * plot_prcurve
    * plot_pulse_comparison
    * plot_trifid_appris
    * plot_transcript_types_histogram
    * plot_validation_curve    
"""
from __future__ import absolute_import, division, print_function

import pandas as pd
import numpy as np
import altair as alt
import matplotlib.pyplot as plt
import shap
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_curve, auc, matthews_corrcoef
from yellowbrick.model_selection import LearningCurve, ValidationCurve
from yellowbrick.style import set_palette


# matplotlib
title_font = {'fontname':'Arial', 'size':16, 'color':'black', 'weight':'normal', 'verticalalignment':'bottom'}
axis_font = {'fontname':'Arial', 'size':14}

def cat_appris_order(df:pd.DataFrame) -> pd.DataFrame:
    """It creates an way of sort the APPRIS histogram.

    Args:
        df (pd.DataFrame): Labels unsorted.

    Returns:
        pd.DataFrame: Labels sorted.
    """    
    df.loc[df['appris'] == 'PRINCIPAL:1', 'appris_order'] = 1
    df.loc[df['appris'] == 'PRINCIPAL:2', 'appris_order'] = 2
    df.loc[df['appris'] == 'PRINCIPAL:3', 'appris_order'] = 3
    df.loc[df['appris'] == 'PRINCIPAL:4', 'appris_order'] = 4
    df.loc[df['appris'] == 'PRINCIPAL:5', 'appris_order'] = 5
    df.loc[df['appris'] == 'ALTERNATIVE:1', 'appris_order'] = 6
    df.loc[df['appris'] == 'ALTERNATIVE:2', 'appris_order'] = 7
    df.loc[df['appris'] == 'MINOR', 'appris_order'] = 8
    return df


def cat_transcript_type(df):
    '''
    It creates a category name for transcript types.
    '''
    df.loc[df['flags'].str.contains('IG'), 'flags_mod'] = 'Immunoglobulin'
    df.loc[df['flags'].str.contains('TR'), 'flags_mod'] = 'T-cell receptor'
    df.loc[df['flags'].str.contains('protein_coding'), 'flags_mod'] = 'Protein coding'
    df.loc[df['flags'].str.contains('nonsense_mediated_decay'), 'flags_mod'] = 'Nonsense mediated decay'
    df.loc[df['flags'].str.contains('non_stop_decay'), 'flags_mod'] = 'Non stop decay'
    df.loc[df['flags'].str.contains('polymorphic_pseudogene'), 'flags_mod'] = 'Polymorphic pseudogenes'
    return df


def config_altair(chart:object, height:int, width:int)->object:
    """It makes a configuration predefined for altair charts.

    Args:
        chart (object): altair chart to input in the function. 
        height (int): figure height.
        width (int): figure width.

    Returns:
        [object]: altair chart updated.
    """    
    return chart.properties(
        height=height,
        width=width,
        title=''
    ).resolve_scale(
        y='independent'
    ).configure_axis(
        labelFontSize=16,
        titleFontSize=18
    )


def create_categories(df_features:list, df_predictions:list, feature:str, cats:list)->list:
    """Create categories for TRIFID vs. Rest visualization

    Args:
        df_features (list): features pandas DataFrame.
        df_predictions (list): predictions pandas DataFrame.
        feature (str): feature to categorize.
        cats (list): categories.

    Returns:
        list: pandas DataFrame categorized.
    """        
    df_predictions_pc = df_predictions[df_predictions['flags'].str.contains('protein_coding')]
    df = pd.merge(
        df_features[['transcript_id', feature]],
        df_predictions_pc[['gene_name', 'transcript_id', 'trifid_score']], 
        how='left', on='transcript_id')
    df_g = df.groupby('gene_name')['trifid_score', feature].max()
    
    df_g.loc[df_g[feature] < cats[0], f'{feature.capitalize()}'] = f'Lower than {cats[0]}'
    df_g.loc[df_g[feature] > cats[-1], f'{feature.capitalize()}'] = f'More than {cats[-1]}'
    
    df_g.loc[df_g[f'{feature.capitalize()}'] == f'Lower than {cats[0]}', 'cat_order'] = 1
    df_g.loc[df_g[f'{feature.capitalize()}'] == f'More than {cats[-1]}', 'cat_order'] = len(cats)+1
    for n, (i,k) in enumerate(zip(cats,cats[1:])):
        df_g.loc[(df_g[feature] >= i) & (df_g[f'{feature}'] <= k), f'{feature.capitalize()}'] = f'{i} to {k}'
        df_g.loc[df_g[f'{feature.capitalize()}'] == f'{i} to {k}', 'cat_order'] = n+2
    df_g['cat_order'] = df_g['cat_order'].astype(int)
    df_g = df_g.sort_values(by='cat_order', ascending=False)
    df_g['feature'] = feature
    return df_g


def explain_prediction(df:pd.DataFrame, model:object, features:list, transcript_id:str)->pd.DataFrame:
    """It explains the prediction of a specific transcript.

    Args:
        df (pd.DataFrame): pandas DataFrame with the features.
        model (object): trained model.
        features (list): features to use in the prediction.
        transcript_id (str): transcript id to explain.

    Returns:
        pd.DataFrame: pandas DataFrame with the explanation.
    """
    trifid_model = model.model
    train_features = model.train_features
    train_target = model.train_target
    trifid_model.fit(train_features, train_target)
    df_i = df.set_index('transcript_id')
    df_features = df_i.iloc[df_i.index.get_level_values('transcript_id') == transcript_id]
    df_features = df_features[features]
    explainer = shap.TreeExplainer(trifid_model)
    shap_values = explainer.shap_values(df_features)
    df_shap_values = pd.DataFrame(
        list(zip(np.abs(shap_values).mean(0)[0]
            )),
        columns=['shap_value'], index=df_features.columns)
    df_shap_values = df_shap_values.sort_values('shap_value', ascending=False)
    
    print('The functional score predicted by the model for {} is {}'.format(
        transcript_id, 
        round(df_i.iloc[df_i.index.get_level_values('transcript_id') == transcript_id]['trifid_score'].values[0], 3)))
    print('The SHAP base value is {}'.format(
        round(explainer.expected_value[1],4)))
    plot = shap.force_plot(explainer.expected_value[1], shap_values[1], df_features, 
                    feature_names=None, out_names=f'value predicted for {transcript_id}', link='identity', 
                    plot_cmap=['#0099bf', '#ee9c00'], matplotlib=True, show=False, figsize=(20, 3), 
                    ordering_keys=True, ordering_keys_time_format=None, text_rotation=None)
                    

def plot_appris_histogram(df:pd.DataFrame)->object:
    '''It creates an histogram with predicted score over APPRIS labels.

    Args:
        df (pd.DataFrame): pandas DataFrame with the features.

    Returns:
        object: altair chart.
    '''
    source = df.sort_values(by='appris_order', ascending=True)
    xcol='trifid_score'
    ycol='appris'
    names_order = ['PRINCIPAL:1', 'PRINCIPAL:2', 'PRINCIPAL:3', 'PRINCIPAL:4',
                   'PRINCIPAL:5', 'ALTERNATIVE:1', 'ALTERNATIVE:2', 'MINOR']
    color_order =['#8fc99e', '#8fc99e', '#8fc99e', '#8fc99e', '#8fc99e', '#F8DFC0', '#F8DFC0', '#b8b8b8']
    opacity=1
    chart = alt.Chart(source[['appris', 'trifid_score']]).mark_bar().encode(
        x=alt.X(xcol,
                bin=alt.Bin(extent=[0, 1],
                            step=0.05),
                axis=alt.Axis(title=f'TRIFID Score', titleFont='Arial', titleFontWeight='normal', titleFontSize=20, 
                            labelFont='Arial', labelFontSize=16)),
        y=alt.Y('count()', axis=alt.Axis(title='', labelFont='Arial', labelFontSize=16)),
        opacity=alt.value(opacity),
        color=alt.Color(ycol,
                        legend=None,
                        title='',
                        scale=alt.Scale(domain=names_order, range=color_order)),
        row=alt.Row(ycol,
                    title='',
                    header=alt.Header(labelAngle=0, labelAlign='left',
                    labelFontSize=16, labelFont='Arial', labelFontWeight='normal'),
                    sort=names_order)
    )
    return config_altair(chart, height=42.5, width=550)


def plot_feature_importances(source, xcol, ycol, facetcol, method, ntop=18):
    '''
    It plots the feature importances of a model spplited by category (facetcol).
    '''
    opacity = 0.8
    source = source.sort_values(by=xcol, ascending=False)[:ntop]
    palette = alt.Scale(
        range=[
            '#BAB0AC', # gray annotation
            '#4E79A7', # blue evolution
            '#59A14F', # green expression
            '#EDC948', # yellow splicing
            '#F28E2B' # orange structure
            ]
        )
    feature_importances = alt.Chart().mark_bar().encode(
        x=alt.X(xcol, axis=alt.Axis(title=f'Top {method} values (% impact on model output)',
                                    titleFont='Arial',
                                    titleFontWeight='normal',
                                    titleFontSize=20)),
        y=alt.Y(ycol, axis=alt.Axis(title='', titleFont='Arial', titleFontWeight='normal', grid=True),
                type='nominal', sort=alt.EncodingSortField(field=xcol, order='descending'),
                scale=alt.Scale()
                ),
        opacity=alt.value(opacity),
        color=alt.Color(facetcol,
                        legend=None,
                        title='',
                        scale=palette)
    )

    chart = alt.layer(feature_importances, data=source).facet(
        row=alt.Row(facetcol,
                    title='',
                    header=alt.Header(labelAngle=0, labelFontSize=16, labelFont='Arial',
                    labelOrient='right', labelPadding=20),
                    sort = alt.EncodingSortField('name',  order='ascending')
                    )
    ).resolve_scale(
        x='shared',
        y='independent',
    ).configure_scale(
        bandPaddingInner=0.2,
    ).configure_facet(
        spacing=10,
    ).configure_view(
        width=750
    ).configure_axisY(
        labelColor='black',
        domainWidth=1.5,
        domainOpacity=1.5,
        domainColor='gray',
        labelFontSize=14,
        labelFont='Arial'
    ).configure_axisX(
        domainWidth=1.5,
        domainOpacity=1.5,
        domainColor='gray',
        labelFontSize=16,
        labelFont='Arial'
    )
    return chart


def plot_learning_curve(model:object, seed:int, train_size:list, scoring_function:object, features:list, target:list, ylabel:str, ylim=list):
    """A learning curve shows the relationship of the training score versus the 
    cross validated test score for an estimator with a varying number of 
    training samples. 

    https://www.scikit-yb.org/en/latest/api/model_selection/learning_curve.html

    Args:
        model (object): scikit-learn model to be evaluated.
        seed (int): random state.
        train_size (list): training instances used to generate the learning curve.
        scoring_function (object): scikit-learn scoring function.
        features (list): pandas DataFrame with features to train the model.
        target (list): pandas DataFrame with labels to train the model.
        ylabel(str): name of the metric.
        ylim(list): lower and upper limit.
    """    
    set_palette("dark")
    fig = plt.figure()                                                                                                                     
    ax = fig.add_subplot(111)
    cv = StratifiedKFold(n_splits=5, shuffle=False, random_state=seed)
    viz = LearningCurve(model, cv=cv, scoring=scoring_function, train_sizes=train_size, 
                        n_jobs=-1, size=(500, 350))
    viz.fit(features, target)    
    viz.ax.set_xlabel("Training samples", fontsize=16)                                                                                                 
    viz.ax.set_ylabel(ylabel, fontsize=16) 
    viz.ax.legend(loc='lower right', fontsize=16, bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    viz.ax.tick_params(labelsize=14)
    viz.ax.grid(True)
    viz.ax.set_ylim(ylim)
    viz.set_title('')
    return fig
    

def plot_prcurve(model:object, X:list, y:list, n_splits:str, seed:int, ax:object, title:str=None):
    """It plots a classic precision recall curve

    Args:
        model (object): scikit-learn model to be evaluated.
        X (list): pandas DataFrame with features to train the model.
        y (list): pandas DataFrame with labels to train the model.
        n_splits (str): CV splits.
        seed (int): random state.
        ax (object): matplotlib axis.
        title (str, optional): plot title. Defaults to None.
    """    
    k_fold = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    y_real = []
    y_pred = []
    y_proba = []
    for i, (train_index, test_index) in enumerate(k_fold.split(X, y)):
        Xtrain, Xtest = X[train_index], X[test_index]
        ytrain, ytest = y[train_index], y[test_index]
        model.fit(Xtrain, ytrain)
        preds = model.predict(Xtest)
        pred_proba = model.predict_proba(Xtest)
        precision, recall, _ = precision_recall_curve(ytest, pred_proba[:,1])
        #lab = 'Fold %d AUPR=%.3f' % (i+1, auc(recall, precision))
        ax.step(recall, precision, lw=1.3, color='#636363', alpha=.3, where='post')
        y_real.append(ytest)
        y_pred.append(preds)
        y_proba.append(pred_proba[:,1])

    y_real = np.concatenate(y_real)
    # y_pred = np.concatenate(y_pred)
    y_proba = np.concatenate(y_proba)
    precision, recall, _ = precision_recall_curve(y_real, y_proba)
    lab = 'Overall AUC-PR=%.3f' % (auc(recall, precision))
    # print(f'MCC = {matthews_corrcoef(y_real, y_pred)}')
    ax.step(recall, precision, label=lab, lw=1.8, color='#de425b', where='post')
    ax.set_xlabel('Recall', family='Arial', fontsize=16)
    ax.set_ylabel('Precision', family='Arial', fontsize=16)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, 1)
    ax.legend(loc='best', fontsize=14)
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)
    ax.grid(False)
    ax.tick_params(labelsize=13)
    if title:
        ax.set_title(title, fontsize=16)
    return ax


def plot_pulse_comparison(df, x_col='trifid_score', x_axis_name='TRIFID score'):
    '''
    It plots a binned 2d heatmap which compares the distribution with PULSE scores.
    '''
    base = alt.Chart(df)

    xscale = alt.Scale(domain=(0, 1))
    yscale = alt.Scale(domain=(0, 1))
    area_args = {'opacity': 0.8,
                 'interpolate': 'step',
                 'color': '#2F4F4F',
                 'limit': 1}
    bins = [0, 0.2, 0.4, 0.6, 0.8, 1]
    heatmap = base.mark_rect().encode(
        alt.X(x_col, bin=alt.Bin(maxbins=10),
              axis=alt.Axis(title=x_axis_name,
                            titleFont='Arial',
                            titleFontWeight='normal',
                            titleFontSize=28,
                            titlePadding=15,
                            labelFontSize=26,
                            labelFont='Arial',
                            labelFontWeight='normal',
                            format='.1f',
                            values=bins)),
        alt.Y('PULSE', bin=alt.Bin(maxbins=10),
              axis=alt.Axis(title='PULSE predicted score',
                            titleFont='Arial',
                            titleFontWeight='normal',
                            titleFontSize=28,
                            titlePadding=15,
                            labelFontSize=26,
                            labelFont='Arial',
                            labelFontWeight='normal',
                            format='.1f',
                            values=bins)),
        alt.Color('count()', scale=alt.Scale(scheme='teals'))
    ).properties(
        width=400,
        height=400
    )

    top_hist = base.mark_bar(**area_args).encode(
        alt.X(x_col,
              bin=alt.Bin(maxbins=20, extent=xscale.domain),
              stack=None,
              title='',
              axis=alt.Axis(values=[])
              ),
        alt.Y('count()', stack=None, title='', axis=alt.Axis(
            labelFontSize=26, labelFont='Arial', labelFontWeight='normal')),
        # alt.Color('transcript_type:N'),
    ).properties(
        height=50,
        width=400
    )

    right_hist = base.mark_bar(**area_args).encode(
        alt.Y('PULSE:Q',
              bin=alt.Bin(maxbins=20, extent=yscale.domain),
              stack=None,
              title='',
              axis=alt.Axis(values=[])
              ),
        alt.X('count()',
              stack=None,
              title='',
              axis=alt.Axis(
            labelFontSize=26, labelFont='Arial', labelFontWeight='normal')
              ),
        # alt.Color('transcript_type:N'),
    ).properties(
        width=50,
        height=400
    )
    chart = (top_hist & (heatmap | right_hist)).configure(
        concat=alt.CompositionConfig(spacing=5),
        countTitle='Number of isoforms'
    ).configure_legend(
        gradientLength=400,
        gradientOpacity=1,
        gradientThickness=20,
        orient='bottom',
        gradientDirection='horizontal',
        titleFontSize=20,
        titleFont='Arial',
        titleFontWeight='normal',
        tickCount=6,
        labelFontSize=22,
        labelFont='Arial',
        labelColor='#000000'
    ).configure_rect(
        opacity=0.8,
        binSpacing=0.5
    )
    return chart


def plot_trifid_appris(df):
    '''
    It plots a paired histogram comparin PRINCIPAL with ALTERNATIVE isoforms.
    '''
    xcol='trifid_score'
    ycol='appris'
    df = df[[xcol, ycol]]

    df.loc[df['appris'].str.contains('PRINCIPAL'), 'appris_label'] = 'PRINCIPAL'
    df.loc[~df['appris'].str.contains('PRINCIPAL'), 'appris_label'] = 'ALTERNATIVE'

    opacity=0.95
    names_order = ['PRINCIPAL', 'ALTERNATIVE']
    color_order = ['#A8C48A', '#EEC99F']

    chart = alt.Chart(df).transform_joinaggregate(
        total='count(*)'
    ).transform_calculate(
        pct='1 / datum.total'
    ).mark_bar().encode(
        x=alt.X(xcol,
                bin=alt.Bin(extent=[0, 1],
                            step=0.1),
                axis=alt.Axis(title=f'TRIFID score (binned)', titleFontSize=18, titleFont='Arial', titleFontWeight='normal',
                             labelFontSize=16)),
        y=alt.Y('sum(pct):Q',
                scale=alt.Scale(domain=(0, 0.45)),
                axis=alt.Axis(title='', labelFontSize=16)),
        opacity=alt.value(opacity),
        color=alt.Color('appris_label',
                        legend=None,
                        title='',
                        scale=alt.Scale(domain=names_order, range=color_order)),
        row=alt.Row('appris_label',
                    title='',
                    header=alt.Header(labelAngle=0, labelFontSize=16, labelAlign='left'),
                    sort='descending'))
    return config_altair(chart, height=100, width=600)


def plot_transcript_types_histogram(df):
    '''
    It creates a faceted histogram per transcript type tags.
    '''

    xcol='trifid_score'
    ycol='flags_mod'
    names_order = ['Protein coding', 'Nonsense mediated decay', 'Non stop decay',
                   'Polymorphic pseudogenes', 'Immunoglobulin', 'T-cell receptor']
    palette='set2'
    opacity=1

    chart = alt.Chart(df[['trifid_score', 'flags_mod']]).mark_bar().encode(
        x=alt.X(xcol,
                bin=alt.Bin(extent=[0, 1],
                            step=0.05),
                axis=alt.Axis(title='TRIFID Score', 
                         titleFont='Arial',
                         titleFontWeight='normal',
                         titleFontSize=22,labelFont='Arial', labelFontSize=16)),
        y=alt.Y('count()', axis=alt.Axis(title='', 
                         titleFont='Arial',
                         titleFontWeight='normal',
                         titleFontSize=18,labelFont='Arial', labelFontSize=16)),
        opacity=alt.value(opacity),
        color=alt.Color(ycol,
                        legend=None,
                        title='',
                        scale=alt.Scale(scheme=palette)),
        row=alt.Row(ycol,
                    title='',
                    header=alt.Header(labelAngle=0, 
                    labelFontSize=18, labelFont='Arial', labelFontWeight='normal', labelAlign='left'),
                    sort=names_order)
    )
    return config_altair(chart, height=45, width=500)


def plot_validation_curve(model:object, seed:int, param:str, param_range:list, scoring_function:object, features:list, target:list):
    """Model validation is used to determine how effective an estimator is on 
    data that it has been trained on as well as how generalizable it is to new 
    input.

    https://www.scikit-yb.org/en/latest/api/model_selection/validation_curve.html

    Args:
        model (object): scikit-learn model to be evaluated.
        seed (int): random state.
        param (str): Paramater to be evaluated.
        param_range (list): Range to evaluate the parameter selected.
        scoring_function (object): scikit-learn scoring function.
        features (list): pandas DataFrame with features to train the model
        target (list): pandas DataFrame with labels to train the model
    """    
    fig = plt.figure()                                                                                                                     
    ax = fig.add_subplot(111)
    set_palette('paired')                                                                                                                    
    cv = StratifiedKFold(n_splits=5, shuffle=False, random_state=seed)
    viz = ValidationCurve(model, param_name=param, param_range=param_range, cv=cv, scoring=scoring_function, n_jobs=-1)
    viz.fit(features, target) 
    # plt.savefig(".", bbox_inches='tight')
    plt.show()
    plt.gcf().clear()