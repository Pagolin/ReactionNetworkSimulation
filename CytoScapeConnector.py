# HTTP Client for Python
import requests
from py2cytoscape.data.cyrest_client import CyRestClient
import igraph
"""
THis would be the way to go without py2cytoscape
# Standard JSON library
import json
# Basic Setup
PORT_NUMBER = 8080
BASE = 'http://localhost:' + str(PORT_NUMBER) + '/v1/'
# Header for posting data to the server as JSON
HEADERS = {'Content-Type': 'application/json'}
# Define dictionary of empty network
empty_network = {
        'data': {
            'name': 'I\'m empty!'
        },
        'elements': {
            'nodes':[],
            'edges':[]
        }
}
res = requests.post(BASE + 'networks?collection=My%20Collection', data=json.dumps(empty_network), headers=HEADERS)
new_network_id = res.json()['networkSUID']
print('Empty network created: SUID = ' + str(new_network_id))
"""
cy = CyRestClient()
def plot_cool_network(nx_graph):
    network =  cy.network.create_from_networkx(nx_graph, collection='Generated by NetworkX')
    print(network.get_id())
