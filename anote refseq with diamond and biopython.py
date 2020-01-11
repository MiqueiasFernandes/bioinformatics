#!/usr/bin/env python
# coding: utf-8



## author: Miquéias Fernandes 01/2020  bio@mikeias.net



# In[1]:


uniprot = [l.strip().split('\t') for l in open('uniprot.tab').readlines()]


# In[2]:


uniprot_vir = {y[0]: y[1:] for y in [x for x in uniprot if 'Viridiplantae' in x[7]]}


# In[3]:


uniprot_oth = {y[0]: y[1:] for y in [x for x in uniprot if not 'Viridiplantae' in x[7]]}


# In[4]:


diamond_on_unip = {x# In[5]:


diamond_on_unip_vir = {z[0]:['SP', uniprot_vir[z[1]]] for z in [x for x in diamond_on_unip.items() if x[1] in uniprot_vir]}


# In[6]:


diamond_on_unip_oth = {z[0]:['SO', uniprot_oth[z[1]]] for z in [x for x in diamond_on_unip.items() if x[1] in uniprot_oth]}


# In[7]:


diamond_on_nr = {x[0]: x[1] for x in [l.strip().split('\t') for l in open('../Galaxy19-[Diamond_on_data_16_and_data_11].tabular').readlines()]}


# In[8]:


diamond_on_nr_filt = {y[0]: y[1] for y in [x for x in diamond_on_nr.items() if not x[0] in diamond_on_unip_vir]}


# In[9]:


nr_ids_anot = list(set(diamond_on_nr_filt.values()))


# In[10]:

[0]: x[1].split('|')[1] for x in [l.strip().split('\t') for l in open('../Galaxy13-[Diamond_on_data_12_and_data_11].tabular').readlines()]}



from Bio import Entrez
Entrez.email = "bio@ufes.br"


# In[11]:




def ncbiRefseq(ids, db="protein", chunk=100, retry=True, startsAt=1):
    res = {}
    chunks = [ids[i:min(len(ids),i+chunk)] for i in range(0,len(ids),chunk)]
    for i in range(startsAt-1, len(chunks)):
        try:
            if len(chunks) > 1:
                print('process chunk %d of %d ...' % (i+1, len(chunks)))
            with Entrez.esummary(db=db, id=','.join(chunks[i]), retmode="xml") as handle:
                for r in Entrez.parse(handle):
                        res.update({r['AccessionVersion']: [r['Id'], r['Title'], str(r['TaxId'])]})
        except:
            if retry:
                print('chunk %d falhou tentando %d restantes ...' % (i+1, len(chunks)-i))
                ncbiRefseq(ids, db=db, chunk=chunk, retry=retry, startsAt=i+2)
    return res


def ncbiTaxonomy(ids, chunk=100):
    res = {}
    for qs in [ids[i:min(len(ids),i+chunk)] for i in range(0,len(ids),chunk)]:
        with Entrez.efetch(db="taxonomy", id=",".join(qs), rettype="lineage", retmode="xml") as handle:
            res.update({r['TaxId']:{'ScientificName':r['ScientificName'],'Lineage': r['Lineage']}
                   for r in Entrez.parse(handle)})
    return res



def anotRefSeq(ids):
    anot = ncbiRefseq(ids, chunk=1000)
    taxs = ncbiTaxonomy(list(set([x[2] for x in anot.values()])))
    return {x: {'Id': anot[x][0], 'Title': anot[x][1], 'Taxonomy': [anot[x][2], taxs[anot[x][2]]]} for x in anot}


# In[12]:


nr_anot = anotRefSeq(nr_ids_anot)


# In[16]:


nr_plants = {z[0]: ['NP', [z[1], nr_anot[z[1]]]] for z in [x for x in diamond_on_nr_filt.items() if 
               x[1] in nr_anot and 'Viridiplantae' in nr_anot[x[1]]['Taxonomy'][1]['Lineage']]}


nr_oth = {z[0]: ['NO', [z[1], nr_anot[z[1]]]] for z in [x for x in diamond_on_nr_filt.items() if 
               x[1] in nr_anot and not 'Viridiplantae' in nr_anot[x[1]]['Taxonomy'][1]['Lineage']]}


# In[17]:


diamond_on_tr = {x[0]: x[1].split('|')[1] for x in [l.strip().split('\t') for l in open('../Galaxy20-[Diamond_on_data_18_and_data_11].tabular').readlines()]}


# In[18]:


for x in [x for x in diamond_on_unip_vir.keys()] + [x for x in nr_plants]:
    if x in diamond_on_tr:
        del diamond_on_tr[x]


# In[19]:


open('tr_ids','w').write('\n'.join(set(diamond_on_tr.values())))


# In[20]:


tr = [l.strip().split('\t') for l in open('tr_anotado.tab').readlines()]


# In[21]:


tr_vir = {y[0]: y[1:] for y in [x for x in tr if len(x) > 7 and 'Viridiplantae' in x[7]]}


# In[22]:


tr_ot = {y[0]: y[1:] for y in [x for x in tr if len(x) <= 7 or not 'Viridiplantae' in x[7]]}


# In[23]:


diamond_on_tr_vir = {z[0]:['TP', tr_vir[z[1]]] for z in [x for x in diamond_on_tr.items() if x[1] in tr_vir]}


# In[24]:


diamond_on_tr_ot = {z[0]:['TO', tr_ot[z[1]]] for z in [x for x in diamond_on_tr.items() if x[1] in tr_ot]}


# In[25]:


resultados = {}
resultados.update(diamond_on_tr_ot)
resultados.update(nr_oth)
resultados.update(diamond_on_unip_oth)
resultados.update(diamond_on_tr_vir)
resultados.update(nr_plants)
resultados.update(diamond_on_unip_vir)



# In[26]:


allP = set([x for x in diamond_on_unip.keys()] + [x for x in diamond_on_tr.keys()]+ [x for x in diamond_on_nr.keys()])


# In[27]:


faltam = [x for x in allP if not x in resultados]


# In[28]:


len(faltam)


# In[29]:


len(resultados)


# In[30]:


res_sw = {x[0]: x[1] for x in [x for x in resultados.items() if x[1][0].startswith("S")]}
res_nr = {x[0]: x[1] for x in [x for x in resultados.items() if x[1][0].startswith("N")]}
res_tr = {x[0]: x[1] for x in [x for x in resultados.items() if x[1][0].startswith("T")]}


# In[31]:


len(res_sw)


# In[32]:


len(res_nr)


# In[33]:


len(res_tr)


# In[ ]:


## https://www.ncbi.nlm.nih.gov/protein/XP_030550848.1
## https://www.uniprot.org/uniprot/Q6XKE6


# In[ ]:





# In[35]:


with open('table_final_anotado.tsv', 'w') as o:
    o.write('guava protein\tmode\tentry\ttaxonomy\ttitle\tgo\n')
    o.write('\n'.join([x[0]+'\t'+x[1][0]+'\t'+ x[1][1][0]+'\t'+x[1][1][4]+'\t'+x[1][1][2]+'\t'+x[1][1][5] for x in res_sw.items()]) + '\n')
    o.write('\n'.join([x[0]+'\t'+x[1][0]+'\t'+ x[1][1][0]+'\t'+x[1][1][4]+'\t'+x[1][1][2]+'\t'+x[1][1][5] for x in res_tr.items()]) + '\n')
    o.write('\n'.join([x[0]+'\t'+x[1][0]+'\t'+ x[1][1][0]+'\t'+x[1][1][1]['Taxonomy'][0]+'\t'+x[1][1][1]['Title']+'\t' for x in res_nr.items()]) + '\n')


# In[44]:


all_tax = list(set([x[1][1]['Taxonomy'][0] for x in res_nr.values()] + [x[1][4] for x in res_sw.values()] + [x[1][4] for x in res_tr.values()]))


# In[45]:


with open('table_final_taxons.tsv', 'w') as o:
    o.write('taxId\tScientificName\tLineage\n')
    o.write('\n'.join([x[0] + '\t' + x[1]['ScientificName'] + '\t' + x[1]['Lineage'] for x in ncbiTaxonomy(all_tax).items()]) + '\n')


# In[ ]:




