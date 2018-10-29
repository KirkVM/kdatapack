from xml.dom import minidom
def getChildElements(node):
	thelist=[]
	for nd in node.childNodes:
		if nd.nodeType==nd.ELEMENT_NODE: thelist.append(nd)
	return thelist

def getNamedNode(theName,nodeList):
	for nd in nodeList:
		if nd.nodeName==theName: return nd
	raise ValueError('no such named node')

def getNamedNodeList(theName,nodeList):
	rv=[]
	for nd in nodeList:
		if nd.nodeName==theName: rv.append(nd)
	return rv
