from interFEBio.XPLT.XPLT import xplt
import numpy as np
import xml.etree.ElementTree as ET
from xml.dom import minidom

file = xplt("jobs/Biaxial.xplt")
file.readAllStates()
disp = file.results.node["displacement"][:, :, :]

xForce = (
    -file.results.face_region["surface reaction force"]
    .region("PrescribedDisplacement2")
    .time(":")
    .comp("x")
)

yForce = (
    file.results.face_region["surface reaction force"]
    .region("PrescribedDisplacement3")
    .time(":")
    .comp("y")
)

meshNodes = file.mesh.nodes.xyz
minX = meshNodes[:, 0].min()
maxX = meshNodes[:, 0].max()
minY = meshNodes[:, 1].min()
maxY = meshNodes[:, 1].max()

ustarx = (meshNodes[:, 0] - minX) / (maxX - minX)
ustary = (-meshNodes[:, 1] + maxY) / (maxY - minY)
ustarz = np.zeros_like(ustarx)

ustarx = np.column_stack((ustarx, ustarz, ustarz)) * 2
ustary = np.column_stack((ustarz, ustary, ustarz)) * 2

root = ET.Element("febio_optimize")
md = ET.SubElement(root, "MeasuredDisplacements")
for t in range(disp.shape[0]):
    time_el = ET.SubElement(md, "time", {"t": str(t)})
    for nodeID in range(disp.shape[1]):
        # print(disp[t, nodeID, :])
        nodeEl = ET.SubElement(time_el, "node", {"id": str(nodeID + 1)})
        nodeEl.text = ",".join("{:.6g}".format(float(x)) for x in disp[t, nodeID])


vdMain = ET.SubElement(root, "VirtualDisplacements")
vd1 = ET.SubElement(vdMain, "virtualdisplacement", {"id": "1"})
time1 = ET.SubElement(vd1, "time", {"t": str(1)})
for nodeID in range(ustarx.shape[0]):
    nodeUstar = ET.SubElement(time1, "node", {"id": str(nodeID + 1)})
    nodeUstar.text = ",".join("{:.6g}".format(float(x)) for x in ustarx[nodeID])


vd2 = ET.SubElement(vdMain, "virtualdisplacement", {"id": "2"})
time1 = ET.SubElement(vd2, "time", {"t": str(1)})
for nodeID in range(ustary.shape[0]):
    nodeUstar = ET.SubElement(time1, "node", {"id": str(nodeID + 1)})
    nodeUstar.text = ",".join("{:.6g}".format(float(x)) for x in ustary[nodeID])


loads = ET.SubElement(root, "MeasuredLoads")
for tforce in range(xForce.shape[0]):
    time = ET.SubElement(loads, "time", {"t": str(tforce + 1)})
    xforceEl = ET.SubElement(time, "surface", {"id": "PrescribedDisplacement2"})
    xforceEl.text = ",".join("{:.6g}".format(float(x)) for x in xForce[tforce])
    yforceEl = ET.SubElement(time, "surface", {"id": "PrescribedDisplacement3"})
    yforceEl.text = ",".join("{:.6g}".format(float(-x)) for x in yForce[tforce])

xml_bytes = ET.tostring(root, encoding="utf-8")
pretty = minidom.parseString(xml_bytes).toprettyxml(indent="  ")

print(pretty)

with open("VFMData.feb", "w", encoding="utf-8") as f:
    f.write(pretty)
