newick = {'label': 'Node1', 'children': [{'label': 'Lilium', 'children': [], 'features': 0.1005098375}, {'label': 'Node2', 'children': [{'label': 'Zea', 'children': [], 'features': 0.0405061803}, {'label': 'Node3', 'children': [{'label': 'Triticum', 'children': [], 'features': 0.0022100542}, {'label': 'Hordeum', 'children': [], 'features': 0.0077439602}], 'features': 0.0226705269}], 'features': 0.0106110864}, {'label': 'Oryza', 'children': [], 'features': 0.0122492691}], 'features': None}

'''print(newick)

print(type(newick))

print(newick['children'])'''

label = newick['label']
children = newick['children']
features = newick['features']

done = False
while not done:
    print(label)
    for child in children:
        print(child['label'])
    print(features)




