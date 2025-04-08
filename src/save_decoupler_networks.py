import decoupler as dc

net = dc.get_dorothea()
net.to_csv("resources/decoupler/dorothea.csv", index=False)

net = dc.get_collectri()
net.to_csv("resources/decoupler/collectri.csv", index=False)
