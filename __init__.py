bl_info = {
    "name": "SwiftBlock",
    "author": "Karl-Johan Nogenmyr",
    "version": (0, 1),
    "blender": (2, 6, 6),
    "api": 44000,
    "location": "Tool Shelf",
    "description": "Writes block geometry as blockMeshDict file",
    "warning": "not much tested yet",
    "wiki_url": "http://openfoamwiki.net/index.php/SwiftBlock",
    "tracker_url": "",
    "support": 'COMMUNITY',
    "category": "OpenFOAM"}

#----------------------------------------------------------
# File scene_props.py
#----------------------------------------------------------
import bpy
from bpy.props import *

def getPolyLines(verts, edges, obj):
    scn = bpy.context.scene
    polyLinesPoints = []
    polyLines = ''
    polyLinesLengths = [[], []]
    
    def isPointOnEdge(point, A, B):
        eps = (((A - B).magnitude - (point-B).magnitude) - (A-point).magnitude)
        return True if (abs(eps) < scn.tol) else False

    nosnap= [False for i in range(len(edges))]
    for eid, e in enumerate(obj.data.edges):
        nosnap[eid] = e.use_edge_sharp

    bpy.ops.wm.context_set_value(data_path="tool_settings.mesh_select_mode", value="(True,False,False)")
    geoobj = bpy.data.objects[scn.geoobjName]
    geo_verts = list(blender_utils.vertices_from_mesh(geoobj))
    geo_edges = list(blender_utils.edges_from_mesh(geoobj))
    geoobj.select = False # avoid deletion

# First go through all vertices in the block structure and find vertices snapped to edges
# When found, add a vertex at that location to the polyLine object by splitting the edge
# Create a new Blender object containing the newly inserted verts. Then use Blender's
# shortest path algo to find polyLines.

    for vid, v in enumerate(verts):
        found = False
        for gvid, gv in enumerate(geo_verts):
            mag = (v-gv).magnitude
            if mag < scn.tol:
                found = True
                break   # We have found a vertex co-located, continue with next block vertex
        if not found:
            for geid, ge in enumerate(geo_edges):
                if (isPointOnEdge(v, geo_verts[ge[0]], geo_verts[ge[1]])):
                    geo_verts.append(v)
                    geo_edges.append([geo_edges[geid][1],len(geo_verts)-1]) # Putting the vert on the edge, by splitting it in two.
                    geo_edges[geid][1] = len(geo_verts)-1
                    break # No more iteration, go to next block vertex

    mesh_data = bpy.data.meshes.new("deleteme")
    mesh_data.from_pydata(geo_verts, geo_edges, [])
    mesh_data.update()
    geoobj = bpy.data.objects.new('deleteme', mesh_data)
    bpy.context.scene.objects.link(geoobj)
    geo_verts = list(blender_utils.vertices_from_mesh(geoobj))
    geo_edges = list(blender_utils.edges_from_mesh(geoobj))
    bpy.context.scene.objects.active=geoobj

# Now start the search over again on the new object with more verts
    snapped_verts = {}
    for vid, v in enumerate(verts):
        for gvid, gv in enumerate(geo_verts):
            mag = (v-gv).magnitude
            if mag < scn.tol:
                snapped_verts[vid] = gvid
                break   # We have found a vertex co-located, continue with next block vertex
    
    bpy.ops.wm.context_set_value(data_path="tool_settings.mesh_select_mode", value="(True,False,False)")
    for edid, ed in enumerate(edges):
        if ed[0] in snapped_verts and ed[1] in snapped_verts and not nosnap[edid]:
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
            polyLineobj = bpy.data.objects['deleteme.001']
            if len(polyLineobj.data.vertices) > 2:
                polyLineverts = list(blender_utils.vertices_from_mesh(polyLineobj))
                polyLineedges = list(blender_utils.edges_from_mesh(polyLineobj))
                for vid, v in enumerate(polyLineverts):
                    mag = (v-verts[ed[0]]).magnitude
                    if mag < scn.tol:
                        startVertex = vid
                        break
                polyLineStr, vectors, length = sortedVertices(polyLineverts,polyLineedges,startVertex)
                polyLinesPoints.append([ed[0],ed[1],vectors])
                polyLinesLengths[0].append([min(ed[0],ed[1]), max(ed[0],ed[1])]) # write out sorted
                polyLinesLengths[1].append(length) 
                polyLine = 'polyLine {} {} ('.format(*ed)
                polyLine += polyLineStr
                polyLine += ')\n'
                polyLines += polyLine

            geoobj.select = False
            polyLineobj.select = True
            bpy.ops.object.delete()
    geoobj.select = True
    bpy.ops.object.delete()
    return polyLines, polyLinesPoints, polyLinesLengths

def sortedVertices(verts,edges,startVert):
    sorted = []
    vectors = []
    sorted.append(startVert)
    vert = startVert
    length = len(edges)+1
    for i in range(len(verts)):
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
    length = 0.
    for vid, v in enumerate(sorted):
        polyLine += '({} {} {})'.format(*verts[v])
        vectors.append(verts[v])
        if vid>=1:
            length += (vectors[vid] - vectors[vid-1]).magnitude
    return polyLine, vectors, length

def patchColor(patch_no):
    color = [(1.0,0.,0.), (0.0,1.,0.),(0.0,0.,1.),(0.707,0.707,0),(0,0.707,0.707),(0.707,0,0.707)]
    return color[patch_no % len(color)]
    
def initProperties():

    bpy.types.Scene.tol = FloatProperty(
        name = "tol", 
        description = "Snapping tolerance",
        default = 1e-6,
        min = 0.)
        
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
        name = "# cells", 
        description = "For forcing the number of cells on edge (0 disables)",
        default = 0,
        min = 0)
 
    bpy.types.Scene.grading = FloatProperty(
        name = "Grading", 
        description = "The size ratio of first and last cell on edge",
        default = 1.0)
 
    bpy.types.Scene.whichCell = EnumProperty(
        items = [('Coarse', 'Coarse', 'Let the coarse cells have the target resolution'), 
                 ('Fine', 'Fine', 'Let the fine cells have the target resolution')
                 ],
        name = "Cell resolution")
        
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
    
    bpy.types.Scene.snapping = EnumProperty(
        items = [('yes', 'Yes', 'The edge gets a polyLine if its vertices are snapped'), 
                 ('no', 'No', 'The edge will be meshed as a straight line')
                 ],
        name = "Edge snapping")

    bpy.types.Scene.removeInternal = BoolProperty(
        name = "Remove internal faces", 
        description = "Should internal faces be removed?",
        default = False)
        
    bpy.types.Scene.createBoundary = BoolProperty(
        name = "Create boundary faces", 
        description = "Should boundary faces be created?",
        default = False)
        
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
            try:
                obj['swiftBlockObj']
                layout.operator("delete.preview")
            except:
                layout.operator("enable.swiftblock")
        else:
            layout = layout.column()
            layout.operator("write.bmdfile")
            layout.operator("create.preview")
            layout.operator("find.broken")
            layout.prop(scn, 'ctmFloat')
            layout.prop(scn, 'resFloat')
            box = layout.box()
            box = box.column()
            box.label(text='Edge settings')
            box.prop(scn, 'setEdges')
            if scn.setEdges:
                try:
                    geoojb = bpy.data.objects[scn.geoobjName]
                    textstr = "Fetching egde's polyLines from " + geoojb.name
                    box.operator("change.geoobj", text=textstr, emboss=False)
#                    box.prop(scn, 'tol') # the tolerance setting do not behave as I expected... do not adjust for now
                    split = box.split()
                    col = split.column()
                    col.operator('nosnap.edge', text='Curved').snapping = True
                    col = split.column()
                    col.operator('nosnap.edge', text='Straight')
                    box.separator()
                except:
                    box.prop(scn, 'geoobjName')

            split = box.split()
            col = split.column()
            col.prop(scn, 'resForce')
            col.prop(scn, 'grading')
            col.operator("flip.edge")
            col = split.column()
            col.operator("set.edgeres")
            col.row().prop(scn,"whichCell", expand=True)
            col.operator("set.grading")
            col.operator("show.edgeinfo", text="Edge settings")

            box = layout.box()
            box = box.column()
            box.label(text='Patch settings')
            box.prop(scn, 'patchName')
            box.prop(scn, 'bcTypeEnum')
            box.operator("set.patchname")
            for m in obj.data.materials:
                try:
                    patchtype = str(' ' + m['patchtype'])
                    split = box.split(percentage=0.2, align=True)
                    col = split.column()
                    col.prop(m, "diffuse_color", text="")
                    col = split.column()
                    col.operator("set.getpatch", text=m.name + patchtype, emboss=False).whichPatch = m.name
                except:
                    pass
            box.operator("repair.faces")

            group = obj.vertex_groups.active
            rows = 2
            if group:
                rows = 4

            layout.label('Block\'s name settings')
            row = layout.row()
            row.template_list("MESH_UL_vgroups", "", obj, "vertex_groups", obj.vertex_groups, "active_index", rows=rows)
            col = row.column(align=True)
            col.operator("object.vertex_group_add", icon='ZOOMIN', text="")
            col.operator("object.vertex_group_remove", icon='ZOOMOUT', text="").all = False
            if group:
                col.separator()
                col.operator("object.vertex_group_move", icon='TRIA_UP', text="").direction = 'UP'
                col.operator("object.vertex_group_move", icon='TRIA_DOWN', text="").direction = 'DOWN'

            if group:
                row = layout.row()
                row.prop(group, "name")

            if obj.vertex_groups and obj.mode == 'EDIT':
                row = layout.row()

                sub = row.row(align=True)
                sub.operator("object.vertex_group_assign", text="Assign").new = False
                sub.operator("object.vertex_group_remove_from", text="Remove")
                sub = row.row(align=True)
                sub.operator("object.vertex_group_select", text="Select")
                sub.operator("object.vertex_group_deselect", text="Deselect")

class OBJECT_OT_edgeInfo(bpy.types.Operator):
    '''Show/edit settings for one selected edge''' 
    bl_idname = "show.edgeinfo"
    bl_label = "Show and edit settings for selected edge"
 
    def execute(self, context):
        bpy.ops.object.mode_set(mode='OBJECT')
        obj = context.active_object
        scn = context.scene
        for e in obj.data.edges:
            if e.select:
                if scn.snapping == 'yes':
                    e.use_edge_sharp = False
                else:
                    e.use_edge_sharp = True
                    
        bpy.ops.set.edgeres()  
        bpy.ops.set.grading()
        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}
 
    def invoke(self, context, event):
        wm = context.window_manager
        obj = context.active_object
        scn = context.scene
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode='EDIT')
        NoSelected = 0
        for e in obj.data.edges:
            if e.select:
                NoSelected += 1
                if e.use_seam:
                    scn.whichCell = 'Fine'
                else:
                    scn.whichCell = 'Coarse'

                if e.bevel_weight >= 0.001:
                    scn.resForce = obj['bevelToResMap'][str(round(e.bevel_weight*100))]
                else:
                    scn.resForce = 0
                    
                if e.crease == 0:
                    scn.grading = 1
                else:
                    scn.grading = obj['creaseToGradMap'][str(round(e.crease*100))]
                
                if e.use_edge_sharp:
                    scn.snapping = 'no'
                else:
                    scn.snapping = 'yes'
        if NoSelected >= 2:
            self.report({'INFO'}, "More than one edge selected!")        
            return{'CANCELLED'}
        elif NoSelected == 0:
            self.report({'INFO'}, "Please select an edge!")        
            return{'CANCELLED'}
        context.window_manager.invoke_props_dialog(self, width=400)
        return {'RUNNING_MODAL'}  
 
    def draw(self, context):
        scn = context.scene
        split = self.layout.split(percentage=0.5)
        col = split.column()
        col.label("Define polyLine:")
        col.label("Grading:")
        col.label("Which cell gets target res:")
        col.label("Forced resolution:")
        col = split.column()

        col.row().prop(scn, "snapping", expand=True)
        col.prop(scn, "grading")
        col.row().prop(scn,"whichCell", expand=True)
        col.prop(scn, "resForce")

class OBJECT_OT_nosnapEdge(bpy.types.Operator):
    '''Force selected edge(s) straight or curved'''
    bl_idname = "nosnap.edge"
    bl_label = "No snap"

    snapping = BoolProperty(default = False)

    def invoke(self, context, event):
        bpy.ops.object.mode_set(mode='OBJECT')
        obj = context.active_object
        NoSelect = 0
        for e in obj.data.edges:
            if e.select:
                NoSelect += 1
                if not self.snapping:
                    e.use_edge_sharp = True
                else:
                    e.use_edge_sharp = False

        if not NoSelect:
            self.report({'INFO'}, "No edge(s) selected!")
            bpy.ops.object.mode_set(mode='EDIT')
            return{'CANCELLED'}

        bpy.ops.object.mode_set(mode='EDIT')
        return {'RUNNING_MODAL'}


class OBJECT_OT_insertSmoother(bpy.types.Operator):
    '''Inserts a smoother'''
    bl_idname = "insert.smoother"
    bl_label = "Insert smoother"

    def execute(self, context):
        try:
            bpy.data.objects[context.scene.geoobjName]
        except:
            self.report({'INFO'}, "Cannot find object for edges!")        
            return{'CANCELLED'}
        import mathutils
        from . import utils
        bpy.ops.object.mode_set(mode='OBJECT')
        scn = context.scene
        obj = context.active_object
        obj.select=False
        geoobj = bpy.data.objects[scn.geoobjName]
        geoobj.hide = False
        centre = mathutils.Vector((0,0,0))
        no_verts = 0
        profile = utils.smootherProfile()
        res = profile.__len__()
        matrix = obj.matrix_world.copy()
        for v in obj.data.vertices:
            if v.select:
                 centre += matrix*v.co 
                 no_verts += 1
        if no_verts == 0:
            self.report({'INFO'}, "Nothing selected!")
            return{'CANCELLED'}

        centre /= no_verts
        if no_verts <= 2:
             centre = mathutils.Vector((0,0,centre[2]))
        for e in obj.data.edges:
            if e.select:
                (v0id, v1id) = e.vertices
                v0 = matrix*obj.data.vertices[v0id].co
                v1 = matrix*obj.data.vertices[v1id].co
                edgevector = v1-v0
                normal = centre-v0
                tang = normal.project(edgevector)
                normal -= tang
                normal.normalize()
                normal *= -edgevector.length
                p = [v0 for i in range(res)]
                e = [[0,1] for i in range(res)]
                for i in range(res):
                    linecoord = float(i)/(res-1)
                    p[i] = (1-linecoord)*v0+linecoord*v1
                    p[i] += 0.05*normal*profile[i]
                for i in range(res-1):
                    e[i] = [i,i+1]
                mesh_data = bpy.data.meshes.new("deleteme")
                mesh_data.from_pydata(p, e, [])
                mesh_data.update()
                addtoobj = bpy.data.objects.new('deleteme', mesh_data)
                bpy.context.scene.objects.link(addtoobj)
                bpy.data.objects['deleteme'].select = True
                geoobj.select = True
                scn.objects.active = geoobj
                bpy.ops.object.join()
                geoobj.select = False

        geoobj.select = True
        scn.objects.active = geoobj
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.remove_doubles(threshold=0.0001, use_unselected=False)
        bpy.ops.object.mode_set(mode='OBJECT')
        geoobj.select = False
        obj.select = True
        scn.objects.active = obj
        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}
    
class OBJECT_OT_flipEdge(bpy.types.Operator):
    '''Flip direction of selected edge(s). This is useful if you want to \
set grading on several edges which are initially misaligned'''
    bl_idname = "flip.edge"
    bl_label = "Flip edge"

    def execute(self, context):
        bpy.ops.object.mode_set(mode='OBJECT')
        obj = context.active_object
        for e in obj.data.edges:
            if e.select:
                (e0, e1) = e.vertices
                e.vertices = (e1, e0)

        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}
    
class OBJECT_OT_deletePreview(bpy.types.Operator):
    '''Delete preview mesh object'''
    bl_idname = "delete.preview"
    bl_label = "Delete preview mesh"

    def execute(self, context):
        bpy.ops.object.mode_set(mode='OBJECT')
        name = ''
        for obj in bpy.data.objects:
            try:
                name = obj['swiftBlockObj']
            except:
                obj.select = False
        bpy.ops.object.delete()
        try:        
            obj = bpy.data.objects[name]
            obj.select = True
            obj.hide = False
            bpy.context.scene.objects.active = obj
        except:
            pass
        return {'FINISHED'}
    
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
        obj.data.use_customdata_edge_crease = True
        obj.data.use_customdata_edge_bevel = True
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
        NoSelected = 0
        for f in obj.data.polygons:
            if f.select:
                NoSelected += 1
        if NoSelected:
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
        else:
            self.report({'INFO'}, "No faces selected!")        
            return{'CANCELLED'}
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
        NoSelect = 0
        for e in obj.data.edges:
            if e.select == True:
                NoSelect += 1
        if not NoSelect:
            self.report({'INFO'}, "No edge(s) selected!")
            bpy.ops.object.mode_set(mode='EDIT')
            return{'CANCELLED'}
            
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
            currentSum += 1 #int(bevelInt)
        if not foundRes:
            bevelToResMap[str(currentSum)] = res  # create a new entry in map
            bevel = currentSum
        for e in obj.data.edges:
            if e.select == True:
                 e.bevel_weight = bevel/100.
        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class OBJECT_OT_SetGrading(bpy.types.Operator):
    '''Set grading on selected edge(s). Use Ctrl-Alt-Space to find out orientation of each edge. \
Cells will be coarser in the edge's z-direction for grading > 1'''
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
        if scn.whichCell == 'Fine':
            use_seam = True
        else:
            use_seam = False
            
        bpy.ops.object.mode_set(mode='OBJECT')
        NoSelect = 0
        for e in obj.data.edges:
            if e.select == True:
                NoSelect += 1
        if not NoSelect:
            self.report({'INFO'}, "No edge(s) selected!")
            bpy.ops.object.mode_set(mode='EDIT')
            return{'CANCELLED'}

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
            currentSum += 1 # int(creaseInt)
        if not foundGrad:
            creaseToGradMap[str(currentSum)] = grad  # create a new entry in map
            crease = currentSum
        for e in obj.data.edges:
            if e.select == True:
                 e.crease = crease/100.
                 e.use_seam = use_seam
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class OBJECT_OT_FindBroken(bpy.types.Operator):
    '''Detect blocks and mark unused edges'''
    bl_idname = "find.broken"
    bl_label = "Diagnose"

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
        verts = list(blender_utils.vertices_from_mesh(obj))
        edges = list(blender_utils.edges_from_mesh(obj))
        refEdges = list(blender_utils.edges_from_mesh(obj))

        log, block_print_out, dependent_edges, face_info, all_edges, faces_as_list_of_nodes = utils.blockFinder(edges, verts, '','', [])
        bpy.ops.wm.context_set_value(data_path="tool_settings.mesh_select_mode", value="(False,True,False)")
        for e in obj.data.edges:
            e.select = True

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
                obj.data.edges[edgeFinder(v0, v1, refEdges)].select = False

        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class OBJECT_OT_RepairFaces(bpy.types.Operator):
    '''Delete internal face and create boundary faces'''
    bl_idname = "repair.faces"
    bl_label = "Repair"

    c = EnumProperty(
        items = [('wall', 'wall', 'Defines the patch as wall'), 
                 ('patch', 'patch', 'Defines the patch as generic patch'),
                 ('empty', 'empty', 'Defines the patch as empty'),
                 ('symmetryPlane', 'symmetryPlane', 'Defines the patch as symmetryPlane'),
                 ],
        name = "Patch type")
        
    def execute(self, context):
        from . import utils
        import imp
        imp.reload(utils)
        from . import blender_utils

        removeInternal = bpy.context.scene.removeInternal
        createBoundary = bpy.context.scene.createBoundary

        if not createBoundary and not removeInternal:
            self.report({'INFO'}, "Repair: Nothing to do!")        
            return{'CANCELLED'}

        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')
        obj = context.active_object
        verts = list(blender_utils.vertices_from_mesh(obj))
        edges = list(blender_utils.edges_from_mesh(obj))

        disabled = []
        bpy.ops.wm.context_set_value(data_path="tool_settings.mesh_select_mode", value="(True,False,False)")
        for group in obj.vertex_groups:
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.select_all(action='DESELECT')
            if group.name == 'disabled':
                bpy.ops.object.vertex_group_set_active(group=group.name)
                bpy.ops.object.vertex_group_select()
                bpy.ops.object.mode_set(mode='OBJECT')
                vlist = []
                for v in obj.data.vertices:
                    if v.select == True:
                        disabled += [v.index]

        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')
        nRemoved, nCreated = utils.repairFaces(edges, verts, disabled, obj, removeInternal, createBoundary)
        self.report({'INFO'}, "Created {} boundary faces and removed {} internal faces".format(nCreated, nRemoved))        
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.invoke_props_dialog(self, width=200)
        return {'RUNNING_MODAL'}  
 
    def draw(self, context):
        scn = context.scene
        obj = context.object
        nPatches = obj.material_slots.__len__()
        self.layout.prop(scn, "removeInternal")
        self.layout.prop(scn, "createBoundary")
        self.layout.label("Assign new faces to patch")
        self.layout.template_list("MATERIAL_UL_matslots", "", obj, "material_slots", obj, "active_material_index", rows=nPatches)


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
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.object.material_slot_select()
        scn.bcTypeEnum = mat['patchtype'] 
        scn.patchName = self.whichPatch
        return {'FINISHED'}

class OBJECT_OT_createPreview(bpy.types.Operator):
    '''Creates a mesh preview as a separate object from selected vertices'''
    bl_idname = "create.preview"
    bl_label = "Preview"
    
    def execute(self, context):
        if context.scene.setEdges:
            try:
                bpy.data.objects[context.scene.geoobjName]
            except:
                self.report({'INFO'}, "Cannot find object for edges!")        
                return{'CANCELLED'}
        import locale
        from . import utils
        import imp, math
        imp.reload(utils)
        from . import blender_utils
        locale.setlocale(locale.LC_ALL, '')
        scn = context.scene
        obj = context.active_object

        bpy.ops.object.mode_set(mode='OBJECT') #update
        bpy.ops.object.mode_set(mode='EDIT')


        verts = list(blender_utils.vertices_from_mesh(obj))
        edges = list(blender_utils.edges_from_mesh(obj))

        toShow = [0 for i in range(len(verts))]  # A block is previewed if all vertices are selected
        allOff = True
        for vid in range(len(verts)):
            if obj.data.vertices[vid].select:
                toShow[vid] = 1
                allOff = False
        if allOff: # All vertices were unselected - proceed with previewing all blocks
            toShow = [1 for i in range(len(verts))]

        disabled = []
        bpy.ops.wm.context_set_value(data_path="tool_settings.mesh_select_mode", value="(True,False,False)")
        for group in obj.vertex_groups:
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.select_all(action='DESELECT')
            if group.name == 'disabled':
                bpy.ops.object.vertex_group_set_active(group=group.name)
                bpy.ops.object.vertex_group_select()
                bpy.ops.object.mode_set(mode='OBJECT')
                vlist = []
                for v in obj.data.vertices:
                    if v.select == True:
                        disabled += [v.index]

        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')
        
        if not allOff:
            for vid, v in enumerate(obj.data.vertices): #restore selection
                if toShow[vid]:
                    v.select = True

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
            gradedEdges.append([[e.vertices[0],e.vertices[1]], grad, e.use_seam])

        bpy.ops.object.mode_set(mode='OBJECT')
        obj.select = False
        if scn.setEdges:
            polyLines, polyLinesPoints, lengths = getPolyLines(verts, edges, obj)
        else:
            polyLinesPoints = []
            lengths = [[], []]
            
        dx0 = scn.resFloat
        effective_lengths = [[], []]
        for e in edges:
            if e in lengths[0]:
                ind = lengths[0].index(e)
                length = lengths[1][ind]
            elif [e[0],e[1]] in lengths[0]:
                ind = lengths[0].index([e[0],e[1]])
                length = lengths[1][ind]
            else:
                length = (verts[e[0]] - verts[e[1]]).magnitude

            if length < 1e-6:
                self.report({'INFO'}, "Zero length edge detected, check block structure!")        
                return{'FINISHED'}
                
            for ge in gradedEdges:
                if ge[0] == e or ge[0] == [e[0],e[1]]:
                    grad = ge[1]
                    fine = ge[2]
                    break
            if fine: #Fine cells should match the target resolution, this means fewer total cells
                if grad > 1.00001:
                    if grad*dx0/length < 1:
                        length *= (math.log(grad)/(math.log(1-dx0/length)-math.log(1-grad*dx0/length)) +1) * dx0/length
                elif grad < 0.99999:
                    if dx0/length < 1:
                        length /= (math.log(grad)/(math.log(1-dx0/length)-math.log(1-grad*dx0/length)) +1) * dx0/length
                else:
                    pass
            else: #Coarse cells should match the target resolution, this means more total cells
                if grad > 1.00001:
                    if dx0/length < 1:
                        length *= (math.log(grad)/(math.log(1-dx0/(length*grad))-math.log(1-dx0/length)) +1) * dx0/length
                elif grad < 0.99999:
                    if dx0/(length*grad) < 1:
                        length /= (math.log(grad)/(math.log(1-dx0/(length*grad))-math.log(1-dx0/length)) +1) * dx0/length
                else:
                    pass
        # Now we have a final effective edge length. Lets store!
            effective_lengths[0].append([min(e[0],e[1]), max(e[0],e[1])]) # write out sorted
            effective_lengths[1].append(length) 

        size, NoCells = utils.preview(edges, verts, toShow, scn.resFloat/scn.ctmFloat, 
                polyLinesPoints, forcedEdges, gradedEdges, effective_lengths, obj.name, obj, disabled)

        cellstr = locale.format("%d", NoCells, grouping=True)
        if not size:
            self.report({'INFO'}, "Preview mesh is empty. Too few vertices selected, or broken block structure!")
        else:
            self.report({'INFO'}, "Cells in mesh: " + cellstr)

        return{'FINISHED'}

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
            default='/opt',
            )
    check_existing = BoolProperty(
            name="Check Existing",
            description="Check and warn on overwriting existing files",
            default=True,
            options={'HIDDEN'},
            )

    logging = BoolProperty(
            name="Log",
            description="Click to enable log files",
            default=False
            )

    def invoke(self, context, event):
        if context.scene.setEdges:
            try:
                bpy.data.objects[context.scene.geoobjName]
            except:
                self.report({'INFO'}, "Cannot find object for edges!")        
                return{'CANCELLED'}
        try:
            self.filepath = context.active_object['path']
        except:
            self.filepath = 'blockMeshDict'
        bpy.context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}
    
    def execute(self, context):
        import locale
        from . import utils
        import imp, math
        imp.reload(utils)
        from . import blender_utils
        locale.setlocale(locale.LC_ALL, '')

        scn = context.scene
        obj = context.active_object
        patchnames = list()
        patchtypes = list()
        patchverts = list()
        patches = list()
        obj['path'] = self.filepath

        bpy.ops.object.mode_set(mode='OBJECT')   # Refresh mesh object
        bpy.ops.object.mode_set(mode='EDIT')

        bpy.ops.wm.context_set_value(data_path="tool_settings.mesh_select_mode", value="(True,False,False)")
        vertexNames = []
        disabled = []
        for group in obj.vertex_groups:
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.select_all(action='DESELECT')
            if not group.name[0] == '_':
                bpy.ops.object.vertex_group_set_active(group=group.name)
                bpy.ops.object.vertex_group_select()
                bpy.ops.object.mode_set(mode='OBJECT')
                vlist = []
                for v in obj.data.vertices:
                    if v.select == True:
                        vlist += [v.index]
                if vlist:
                    if group.name == 'disabled':
                        disabled = vlist
                    else:
                        vertexNames.append([group.name, vlist])

        bpy.ops.object.mode_set(mode='EDIT')
        for mid, m in enumerate(obj.data.materials):
            bpy.ops.mesh.select_all(action='DESELECT')
            obj.active_material_index = mid
            bpy.ops.object.material_slot_select()
            bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.mode_set(mode='EDIT')
            faces = obj.data.polygons
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
            gradedEdges.append([[e.vertices[0],e.vertices[1]], grad, e.use_seam])

        bpy.ops.object.mode_set(mode='OBJECT')
        obj.select = False
        if scn.setEdges:
            polyLines, polyLinesPoints, lengths = getPolyLines(verts, edges, obj)
        else:
            polyLines = ''
            lengths = [[], []]

        dx0 = scn.resFloat
        effective_lengths = [[], []]
        for e in edges:
            if e in lengths[0]:
                ind = lengths[0].index(e)
                length = lengths[1][ind]
            elif [e[0],e[1]] in lengths[0]:
                ind = lengths[0].index([e[0],e[1]])
                length = lengths[1][ind]
            else:
                length = (verts[e[0]] - verts[e[1]]).magnitude

            if length < 1e-6:
                self.report({'INFO'}, "Zero length edge detected, check block structure!")        
                return{'FINISHED'}
                
            for ge in gradedEdges:
                if ge[0] == e or ge[0] == [e[0],e[1]]:
                    grad = ge[1]
                    fine = ge[2]
                    break
            if fine: #Fine cells should match the target resolution, this means fewer total cells
                if grad > 1.00001:
                    if grad*dx0/length < 1:
                        length *= (math.log(grad)/(math.log(1-dx0/length)-math.log(1-grad*dx0/length)) +1) * dx0/length
                elif grad < 0.99999:
                    if dx0/length < 1:
                        length /= (math.log(grad)/(math.log(1-dx0/length)-math.log(1-grad*dx0/length)) +1) * dx0/length
                else:
                    pass
            else: #Coarse cells should match the target resolution, this means more total cells
                if grad > 1.00001:
                    if dx0/length < 1:
                        length *= (math.log(grad)/(math.log(1-dx0/(length*grad))-math.log(1-dx0/length)) +1) * dx0/length
                elif grad < 0.99999:
                    if dx0/(length*grad) < 1:
                        length /= (math.log(grad)/(math.log(1-dx0/(length*grad))-math.log(1-dx0/length)) +1) * dx0/length
                else:
                    pass
        # Now we have a final effective edge length. Lets store!
            effective_lengths[0].append([min(e[0],e[1]), max(e[0],e[1])]) # write out sorted
            effective_lengths[1].append(length) 

        obj.select = True
        bpy.context.scene.objects.active = obj

        NoCells = utils.write(self.filepath, edges, verts, 
            scn.resFloat/scn.ctmFloat, scn.ctmFloat, 
            patches, polyLines, forcedEdges, gradedEdges, 
            effective_lengths, vertexNames, disabled, self.logging)

        cellstr = locale.format("%d", NoCells, grouping=True)
        self.report({'INFO'}, "Cells in mesh: " + cellstr)

        return{'FINISHED'}


initProperties()

def register():
    bpy.utils.register_module(__name__)

def unregister():
    bpy.utils.unregister_module(__name__)

