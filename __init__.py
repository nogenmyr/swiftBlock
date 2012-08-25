bl_info = {
    "name": "SwiftBlock",
    "author": "Karl-Johan Nogenmyr",
    "version": (0, 1),
    "blender": (2, 6, 2),  # PLEASE NOTE: There is a bug in in version 2.63 which affects the shortest path code! If you really want 2.63, find a daily build version.
    "api": 35622,
    "location": "Tool Shelf",
    "description": "Writes block geometry as blockMeshDict file",
    "warning": "not much tested yet",
    "wiki_url": "None",
    "tracker_url": "",
    "support": 'COMMUNITY',
    "category": "OpenFOAM"}

#----------------------------------------------------------
# File scene_props.py
#----------------------------------------------------------
import bpy
from bpy.props import *

def sortedVertices(verts,edges,startVert):
    sorted = []
    sorted.append(startVert)
    vert = startVert
    length = len(edges)+1
    while len(sorted) < length:
        for eid, e in enumerate(edges):
            if vert in e:
                if e[0] == vert:
                    sorted.append(e[1])
                else:
                    sorted.append(e[0])
                edges.pop(eid)
                vert = sorted[-1]
                break
    polyLine = ''
    for v in sorted:
        polyLine += '({} {} {})'.format(*verts[v])
    return polyLine

def patchColor(patch_no):
    color = [(1.0,0.,0.), (0.0,1.,0.),(0.0,0.,1.),(0.707,0.707,0),(0,0.707,0.707),(0.707,0,0.707)]
    return color[patch_no % len(color)]
    
def initProperties():

    bpy.types.Scene.ctmFloat = FloatProperty(
        name = "convertToMeters", 
        description = "Conversion factor: Blender coords to meter",
        default = 1.0,
        min = 0.)
        
    bpy.types.Scene.resFloat = FloatProperty(
        name = "Resolution", 
        description = "The average spatial resolution of generated mesh in meter",
        default = 1.0,
        min = 0)
 
    bpy.types.Scene.resForce = IntProperty(
        name = "Resolution", 
        description = "For forcing the number of cells on edge (0 disables)",
        default = 0,
        min = 0)
 
    bpy.types.Scene.grading = FloatProperty(
        name = "Grading", 
        description = "The size ratio of first and last cell on edge",
        default = 1.0)
 
    bpy.types.Scene.setEdges = BoolProperty(
        name = "Set edges", 
        description = "Should edges be fetched from another object?",
        default = False)
        
    bpy.types.Scene.geoobjName = StringProperty(
        name = "Object", 
        description = "Name of object to get edges from (this box disappears when object is found)",
        default = '')

    bpy.types.Scene.bcTypeEnum = EnumProperty(
        items = [('wall', 'wall', 'Defines the patch as wall'), 
                 ('patch', 'patch', 'Defines the patch as generic patch'),
                 ('empty', 'empty', 'Defines the patch as empty'),
                 ('symmetryPlane', 'symmetryPlane', 'Defines the patch as symmetryPlane'),
                 ],
        name = "Patch type")
    
    bpy.types.Scene.patchName = StringProperty(
        name = "Patch name",
        description = "Specify name of patch (max 31 chars)",
        default = "defaultName")
    
    return

#
#    Menu in UI region
#
class UIPanel(bpy.types.Panel):
    bl_label = "SwiftBlock settings"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "object"
    
    def draw(self, context):
        layout = self.layout
        scn = context.scene
        obj = context.active_object
        settings = context.tool_settings
        try:
            obj['swiftblock']
        except:
            layout.operator("enable.swiftblock")
        else:
            layout.operator("write.bmdfile")
            layout.operator("find.broken")
            layout.prop(scn, 'ctmFloat')
            layout.prop(scn, 'resFloat')
            box = layout.box()
            box.prop(scn, 'setEdges')
            if scn.setEdges:
                try:
                    geoojb = bpy.data.objects[scn.geoobjName]
                    textstr = "Fetching egde's polyLines from " + geoojb.name
                    box.operator("change.geoobj", text=textstr, emboss=False)
                except:
                    box.prop(scn, 'geoobjName')

            split = box.split()
            col = split.column()
            col.prop(scn, 'resForce')
            col.operator("set.edgeres")
            col = split.column()
            col.prop(scn, 'grading')
            col.operator("set.grading")
            box = layout.box()
            box.label(text='Patch settings')
            box.prop(scn, 'patchName')
            box.prop(scn, 'bcTypeEnum')
            box.operator("set.patchname")
            for m in obj.data.materials:
                try:
                    patchtype = str(' ' + m['patchtype'])
                    split = box.split(percentage=0.2, align=False)
                    col = split.column()
                    col.prop(m, "diffuse_color", text="")
                    col = split.column()
                    col.operator("set.getpatch", text=m.name + patchtype, emboss=False).whichPatch = m.name
                except:
                    pass

class OBJECT_OT_ChangeGeoObj(bpy.types.Operator):
    '''Click to change object'''
    bl_idname = "change.geoobj"
    bl_label = "Change"

    def execute(self, context):
        context.scene.geoobjName = ''
        return {'FINISHED'}
    
class OBJECT_OT_Enable(bpy.types.Operator):
    '''Enables SwiftBlock for the active object'''
    bl_idname = "enable.swiftblock"
    bl_label = "Enable SwiftBlock"
    
    def execute(self, context):
        obj = context.active_object
        obj['swiftblock'] = True

        bpy.context.tool_settings.use_mesh_automerge = True

        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.material_slot_remove()
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')
        try:
            mat = bpy.data.materials['defaultName']
            patchindex = list(obj.data.materials).index(mat)
            obj.active_material_index = patchindex
        except: 
            mat = bpy.data.materials.new('defaultName')
            mat.diffuse_color = (0.5,0.5,0.5)
            bpy.ops.object.material_slot_add() 
            obj.material_slots[-1].material = mat
        mat['patchtype'] = 'wall'
        bpy.ops.object.editmode_toggle()  
        bpy.ops.object.material_slot_assign()
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.editmode_toggle()  
        return{'FINISHED'}

class OBJECT_OT_SetPatchName(bpy.types.Operator):
    '''Set the given name to the selected faces'''
    bl_idname = "set.patchname"
    bl_label = "Set name"
    
    def execute(self, context):
        scn = context.scene
        obj = context.active_object
        bpy.ops.object.mode_set(mode='OBJECT')
        namestr = scn.patchName
        namestr = namestr.strip()
        namestr = namestr.replace(' ', '_')
        try:
            mat = bpy.data.materials[namestr]
            patchindex = list(obj.data.materials).index(mat)
            obj.active_material_index = patchindex
        except: # add a new patchname (as a blender material, as such face props are conserved during mesh mods)
            mat = bpy.data.materials.new(namestr)
            mat.diffuse_color = patchColor(len(obj.data.materials)-1)
            bpy.ops.object.material_slot_add() 
            obj.material_slots[-1].material = mat
        mat['patchtype'] = scn.bcTypeEnum
        bpy.ops.object.editmode_toggle()  
        bpy.ops.object.material_slot_assign()
        return {'FINISHED'}

class OBJECT_OT_SetEdgeRes(bpy.types.Operator):
    '''Force a resolution on selected edge(s)'''
    bl_idname = "set.edgeres"
    bl_label = "Force resolution"
# This very messy way to keep track of resolution is needed as we have to store the info
# in edges native properties. Here bevel_weight is used which is a float in range [0,1]
# The float can store approx. 100 different values. By mult. by 100, int in range [0,100]
# is achieved. This int is mapped to user-set resolution with the 'bevelToResMap'

    def execute(self, context):
        scn = context.scene
        obj = context.active_object
        res = scn.resForce
        bpy.ops.object.mode_set(mode='OBJECT')
        try:
            obj['bevelToResMap']
        except:
            obj['bevelToResMap']= {}
        if res == 0:
            for e in obj.data.edges:
                if e.select == True:
                     e.bevel_weight = 0
            bpy.ops.object.mode_set(mode='EDIT')
            return {'FINISHED'}
        bevelToResMap = obj['bevelToResMap']
        foundRes = False
        currentSum = 1
        for bevelInt in bevelToResMap:
            if bevelToResMap[bevelInt] == res:  # previously used resolution - reuse!
                bevel = int(bevelInt)
                foundRes = True
            currentSum += int(bevelInt)
        if not foundRes:
            bevelToResMap[str(currentSum)] = res  # create a new entry in map
            bevel = currentSum
        for e in obj.data.edges:
            if e.select == True:
                 e.bevel_weight = bevel/100.
        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class OBJECT_OT_SetGrading(bpy.types.Operator):
    '''Set grading on selected edge (Ctrl-Alt-Space to find out orientation)'''
    bl_idname = "set.grading"
    bl_label = "Set grading"
# This very messy way to keep track of grading is needed as we have to store the info
# in edges native properties. Here crease is used which is a float in range [0,1]
# The float can store approx. 100 different values. By mult. by 100, int in range [0,100]
# is achieved. This int is mapped to user-set grading with the 'creaseToGradMap'

    def execute(self, context):
        scn = context.scene
        obj = context.active_object
        grad = scn.grading
        bpy.ops.object.mode_set(mode='OBJECT')
        try:
            obj['creaseToGradMap']
        except:
            obj['creaseToGradMap']= {}
        obj.data.show_edge_crease = False # could be set true to show which edges have grading
        if grad == 1:
            for e in obj.data.edges:
                if e.select == True:
                     e.crease = 0
            bpy.ops.object.mode_set(mode='EDIT')
            return {'FINISHED'}

        creaseToGradMap = obj['creaseToGradMap']
        foundGrad = False
        currentSum = 1
        for creaseInt in creaseToGradMap:
            if creaseToGradMap[creaseInt] == grad:  # previously used grading - reuse!
                crease = int(creaseInt)
                foundGrad = True
            currentSum += int(creaseInt)
        if not foundGrad:
            creaseToGradMap[str(currentSum)] = grad  # create a new entry in map
            crease = currentSum
        for e in obj.data.edges:
            if e.select == True:
                 e.crease = crease/100.
        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class OBJECT_OT_FindBroken(bpy.types.Operator):
    '''Detect blocks and mark unused edges'''
    bl_idname = "find.broken"
    bl_label = "Diagnose"

    if "bpy" in locals():
        import imp
        if "utils" in locals():
            imp.reload(utils)
        if "blender_utils" in locals():
            imp.reload(blender_utils)

    def execute(self, context):
        from . import utils
        import imp
        imp.reload(utils)
        from . import blender_utils
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')
        obj = context.active_object
        obj.data.show_edge_sharp = True
        verts = list(blender_utils.vertices_from_mesh(obj))
        edges = list(blender_utils.edges_from_mesh(obj))
        refEdges = list(blender_utils.edges_from_mesh(obj))

        log, block_print_out, dependent_edges, face_info, all_edges, offences, faces_as_list_of_nodes = utils.blockFinder(edges, verts, '','')
        bpy.ops.wm.context_set_value(data_path="tool_settings.mesh_select_mode", value="(False,True,False)")
        for e in obj.data.edges:
            e.use_edge_sharp = True

        def edgeFinder(v0, v1, edgeList):
            if [v0, v1] in edgeList:
                return edgeList.index([v0, v1])
            if [v1, v0] in edgeList:
                return edgeList.index([v1, v0])
            return -1

        edgeOrder = [[0,1], [1,2], [2,3], [0,3], [4,5], [5,6], [6,7], [4,7], [0,4], [1,5], [2,6], [3,7]]
        for vl in block_print_out:
            for e in edgeOrder:
                v0 = vl[e[0]]
                v1 = vl[e[1]]
                obj.data.edges[edgeFinder(v0, v1, refEdges)].use_edge_sharp = False

        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class OBJECT_OT_GetPatch(bpy.types.Operator):
    '''Click to select faces belonging to this patch'''
    bl_idname = "set.getpatch"
    bl_label = "Get patch"
    
    whichPatch = StringProperty()

    def execute(self, context):
        scn = context.scene
        obj = context.active_object
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.wm.context_set_value(data_path="tool_settings.mesh_select_mode", value="(False,False,True)")
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')
        mat = bpy.data.materials[self.whichPatch]
        patchindex = list(obj.data.materials).index(mat)
        obj.active_material_index = patchindex
        bpy.ops.object.editmode_toggle()  
        bpy.ops.object.material_slot_select()
        scn.bcTypeEnum = mat['patchtype'] 
        scn.patchName = self.whichPatch
        return {'FINISHED'}

class OBJECT_OT_writeBMD(bpy.types.Operator):
    '''Writes out a blockMeshDict file for the currently selected object'''
    bl_idname = "write.bmdfile"
    bl_label = "Write"
    
    if "bpy" in locals():
        import imp
        if "utils" in locals():
            imp.reload(utils)
        if "blender_utils" in locals():
            imp.reload(blender_utils)
    
    filepath = StringProperty(
            name="File Path",
            description="Filepath used for exporting the file",
            maxlen=1024,
            subtype='FILE_PATH',
            )
    check_existing = BoolProperty(
            name="Check Existing",
            description="Check and warn on overwriting existing files",
            default=True,
            options={'HIDDEN'},
            )

    def invoke(self, context, event):
        if context.scene.setEdges:
            try:
                bpy.data.objects[context.scene.geoobjName]
            except:
                self.report({'INFO'}, "Cannot find object for edges!")        
                return{'CANCELLED'}
        self.filepath = 'blockMeshDict'
        bpy.context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}
    
    def execute(self, context):
        from . import utils
        import imp
        imp.reload(utils)
        from . import blender_utils
        scn = context.scene
        obj = context.active_object
        patchnames = list()
        patchtypes = list()
        patchverts = list()
        patches = list()

        bpy.ops.object.mode_set(mode='OBJECT')   # Refresh mesh object
        bpy.ops.object.mode_set(mode='EDIT')
        for mid, m in enumerate(obj.data.materials):
            bpy.ops.mesh.select_all(action='DESELECT')
            obj.active_material_index = mid
            bpy.ops.object.material_slot_select()
            bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.mode_set(mode='EDIT')
            try: 
                faces = obj.data.polygons  # Blender 2.63+
            except:
                faces = obj.data.faces  # Blender 2.62-
            for f in faces:
                if f.select and f.material_index == mid:
                    if m.name in patchnames:
                        ind = patchnames.index(m.name)
                        patchverts[ind].append(list(f.vertices))
                    else:
                        patchnames.append(m.name)
                        patchtypes.append(m['patchtype'])
                        patchverts.append([list(f.vertices)])
                        
        for ind,pt in enumerate(patchtypes):
            patches.append([pt])
            patches[ind].append(patchnames[ind])
            patches[ind].append(patchverts[ind])

        verts = list(blender_utils.vertices_from_mesh(obj))
        edges = list(blender_utils.edges_from_mesh(obj))

        forcedEdges = []
        for e in obj.data.edges:
            if e.bevel_weight >= 0.001:
                N = obj['bevelToResMap'][str(round(e.bevel_weight*100))]
                forcedEdges.append([[e.vertices[0],e.vertices[1]], N])

        gradedEdges = []
        for e in obj.data.edges:
            if e.crease == 0:
                grad = 1
            else:
                grad = obj['creaseToGradMap'][str(round(e.crease*100))]
            gradedEdges.append([[e.vertices[0],e.vertices[1]], grad])

        bpy.ops.object.mode_set(mode='OBJECT')
        polyLines = ''
        
        if scn.setEdges:
            bpy.ops.wm.context_set_value(data_path="tool_settings.mesh_select_mode", value="(True,False,False)")
            geoobj = bpy.data.objects[scn.geoobjName]
            geo_verts = list(blender_utils.vertices_from_mesh(geoobj))
            geo_edges = list(blender_utils.edges_from_mesh(geoobj))
            bpy.context.scene.objects.active=bpy.data.objects[scn.geoobjName]
            geoobj = bpy.data.objects[scn.geoobjName]
            snapped_verts = {}
            for vid, v in enumerate(verts):
                for gv in geo_verts:
                    mag = (v-gv).magnitude
                    if mag < 1e-6:
                        snapped_verts[vid] = geo_verts.index(gv)
                        break
            for ed in edges:
                if ed[0] in snapped_verts and ed[1] in snapped_verts:
                    geoobj.hide = False
                    bpy.ops.object.mode_set(mode='EDIT')
                    bpy.ops.mesh.select_all(action='DESELECT')
                    bpy.ops.object.mode_set(mode='OBJECT')
                    geoobj.data.vertices[snapped_verts[ed[0]]].select = True
                    geoobj.data.vertices[snapped_verts[ed[1]]].select = True
                    bpy.ops.object.mode_set(mode='EDIT')
                    bpy.ops.mesh.select_vertex_path(type='EDGE_LENGTH')

                    bpy.ops.object.mode_set(mode='OBJECT')
                    bpy.ops.object.mode_set(mode='EDIT')
                    bpy.ops.mesh.duplicate()
                    bpy.ops.object.mode_set(mode='OBJECT')
                    bpy.ops.object.mode_set(mode='EDIT')
                    bpy.ops.mesh.separate(type='SELECTED')
                    bpy.ops.object.mode_set(mode='OBJECT')
                    polyLineobj = bpy.data.objects[scn.geoobjName+'.001']
                    polyLineverts = list(blender_utils.vertices_from_mesh(polyLineobj))
                    polyLineedges = list(blender_utils.edges_from_mesh(polyLineobj))
                    for vid, v in enumerate(polyLineverts):
                        mag = (v-verts[ed[0]]).magnitude
                        if mag < 1e-6:
                            startVertex = vid
                            break
                    polyLine = 'polyLine {} {} ('.format(*ed)
                    polyLine += sortedVertices(polyLineverts,polyLineedges,startVertex)
                    polyLine += ')\n'
                    print(polyLine)
                    polyLines += polyLine
                    geoobj.select = False
                    obj.select = False
                    polyLineobj.select = True
                    bpy.ops.object.delete()
        obj.select = True
        bpy.context.scene.objects.active = obj

        utils.write(self.filepath, edges, verts, scn.resFloat/scn.ctmFloat, scn.ctmFloat, patches, polyLines, forcedEdges, gradedEdges)        

        return{'FINISHED'}

initProperties()

def register():
    bpy.utils.register_module(__name__)

def unregister():
    bpy.utils.unregister_module(__name__)

