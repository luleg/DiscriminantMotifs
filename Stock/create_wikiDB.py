import os

larticles = os.listdir('MatNetworks')
larticles = list(filter(lambda a:a[:4]=='wiki',larticles))
#larticles = list(filter(lambda a:'cancel' in a,larticles))

lquali = list(filter(lambda a:'qualite' in a,larticles))
lnonneutre = list(filter(lambda a: 'non-neutre' in a,larticles))
lguerreedition = list(filter(lambda a:'guerre-edition' in a,larticles))

lquali=[art.split('.')[0] for art in lquali]
lnonneutre=[art.split('.')[0] for art in lnonneutre]
lguerreedition=[art.split('.')[0] for art in lguerreedition]

with open('liste_wiki_qualite.txt','w') as f:
    f.write('\n'.join(lquali))


with open('liste_wiki_non-neutre.txt','w') as f:
    f.write('\n'.join(lnonneutre))

with open('liste_wiki_guerre-edition.txt','w') as f:
    f.write('\n'.join(lguerreedition))
