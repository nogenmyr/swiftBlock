    #obj = bpy.data.objects[0]
    #obj.data.vertices.add(1)
    #obj.data.vertices[-1].co
    #obj.data.edges.add(1)
    #obj.data.edges[-1].vertices[0] = Int
    #obj.data.edges[-1].vertices[0] = Int
import mathutils

class lineEdge:
# Rip of OpenFOAM's lineEdge class
    def __init__(self, points):
        self.points_ = points

    def position(self, lambd):
        if (lambd < 0 or lambd > 1):
            print(error)
        return self.points_[0] + lambd * (self.points_[-1] - self.points_[0])

    def length(self):
        return (points_[-1] - points_[0]).magnitude
        
        
class polyLine:
# Rip of OpenFOAM's polyLine class
    def __init__(self, points):

        self.points_ = points  # list of mathtils.vectors... MAKE Sore that points given to init is this type
        self.param_ = [0.0]*len(self.points_) # list of doobles

        if (len(self.param_)):
            for i in range(len(self.param_)-1): #i+1 ??
                self.param_[i+1] = self.param_[i] + (self.points_[i+1] - self.points_[i]).magnitude

            # normalize on the interval 0-1
            self.lineLength_ = self.param_[-1]
            for p in self.param_:
                p /= self.lineLength_
            self.param_[-1] = 1.0
        else:
            self.lineLength_ = 0.0

    def points(self):
        return self.points_

    def nSegments(self):
        return len(self.points_)-1

    def localParameter(self, lambd):
        # check endpoints
        if (lambd < 1.e-10):
            lambd = 0
            return 0
        elif (lambd > 1 - 1.e-10):
            lambd = 1
            return nSegments()

        # search table of cumulative distances to find which line-segment
        # we are on. Check the upper bound.

        segmentI = 1
        while (self.param_[segmentI] < lambd):
            segmentI += 1
        segmentI -= 1   # we want the corresponding lower bound

        # the local parameter [0-1] on this line segment
        lambd = (
            ( lambd - self.param_[segmentI] )
          / ( self.param_[segmentI+1] - self.param_[segmentI] ))

        return segmentI

    def position(self, mu):
        # check endpoints
        if (mu < 1.e-10):
            return self.points_[0]
        elif (mu > 1 - 1.e-10):
            return self.points_[-1]

        lambd = mu
        segment = self.localParameter(lambd)
        return self.position_segment(segment, lambd)

    def position_segment(self, segment, mu):
        # out-of-bounds
        if (segment < 0):
            return self.points_.first()
        elif (segment > self.nSegments()):
            return self.points_.last()

        p0 = self.points()[segment]
        p1 = self.points()[segment+1]

        # special cases - no calculation needed
        if (mu <= 0.0):
            return p0
        elif (mu >= 1.0):
            return p1
        else:
            # linear interpolation
            return self.points_[segment] + mu * (p1 - p0)

    def length(self):
        return self.lineLength_

            
class lineDivide:
# Rip of OpenFOAM's lineDivide class
    def __init__(self, cedge, ndiv, xratio): #cedge can be a polyLine or straightLine, xratio is from calcGexp
        self.points_ = [mathutils.Vector((0,0,0))]*(ndiv+1)  #vectorlist
        self.divisions_ = [0.]*(ndiv+1)
        self.divisions_[ndiv] = 1.0;

        # calculate the spacing
        if (xratio == 1.0):
            for i in range(ndiv-1): #i = 1 ... ndiv-1
                self.divisions_[i+1] = float(i+1)/ndiv;
        else:
            for i in range(ndiv-1):# i = 1 ... ndiv-1
                self.divisions_[i+1] = (1.0 - pow(xratio, i+1))/(1.0 - pow(xratio, ndiv));

        # calculate the points
        for i in range(ndiv+1):
            self.points_[i] =  cedge.position(self.divisions_[i]);

    def points(self):
        return self.points_;

    def lambdaDivisions(self): 
        return self.divisions_;


def buildPreviewMesh(obj, corners, vertices_coord, ni, nj, nk, polyLinesPoints, gradList, vertexNo):
#    from . import polyLine, lineEdge, lineDivide
    def vtxLabel(i, j, k, ni, nj):
        return ( i + j * (ni + 1) + k * (ni + 1) * (nj + 1) )

    def calcGexp(expRatio, dim): #dim = ni. nj. nk
        return pow(expRatio, 1.0/(dim - 1)) if (dim>1) else 0.0

    obj.data.vertices.add( ( (ni+1) * (nj+1) ) * 2  + (ni*2+nj*2) * (nk-1))

    p000 = vertices_coord[corners[0]] 
    p100 = vertices_coord[corners[1]]
    p110 = vertices_coord[corners[2]]
    p010 = vertices_coord[corners[3]]

    p001 = vertices_coord[corners[4]]
    p101 = vertices_coord[corners[5]]
    p111 = vertices_coord[corners[6]]
    p011 = vertices_coord[corners[7]]

    p = [[]*3]*12 # init
    w = [[1]*3]*12 # init

    edgeOrder = [[0,1], [3,2], [7,6], [4,5], [0,3], [1,2], [5,6], [4,7], [0,4], [1,5], [2,6], [3,7] ]
    n = [ni]*4 + [nj]*4 + [nk]*4
    for plp in polyLinesPoints:
        for edgeID in range(12):
            if corners[edgeOrder[edgeID][0]] in plp and corners[edgeOrder[edgeID][1]] in plp:
                polyLineObj = polyLine(plp[2])
                lineDivObj = lineDivide(polyLineObj, n[edgeID], calcGexp(gradList[edgeID], n[edgeID]))
                p[edgeID] = lineDivObj.points()
                w[edgeID] = lineDivObj.lambdaDivisions()

    for edgeID in range(12):
        if not p[edgeID]: # was not set above = is a straight line
            p0 = vertices_coord[corners[edgeOrder[edgeID][0]]]
            p1 = vertices_coord[corners[edgeOrder[edgeID][1]]]
            lineObj = lineEdge([p0,p1])
            n[edgeID]
            lineDivObj = lineDivide(lineObj, n[edgeID], calcGexp(gradList[edgeID], n[edgeID]))
            p[edgeID] = lineDivObj.points()
            w[edgeID] = lineDivObj.lambdaDivisions()

    # generate vertices
    #

    for k in range(nk+1):
        for j in range(nj+1):
            for i in range(ni+1):
                if k==0 or k==nk or j==0 or j==nj or i==0 or i==ni:
#                    vertexNo = vtxLabel(i, j, k, ni, nj) # This is for full 3D mesh, as in blockMesh

                    # points on edges
                    edgex1 = p000 + (p100 - p000)*w[0][i] #   0 1
                    edgex2 = p010 + (p110 - p010)*w[1][i] #    3 2
                    edgex3 = p011 + (p111 - p011)*w[2][i] #   7 6 
                    edgex4 = p001 + (p101 - p001)*w[3][i] #    4 5

                    edgey1 = p000 + (p010 - p000)*w[4][j] #    0 3
                    edgey2 = p100 + (p110 - p100)*w[5][j] #   1 2
                    edgey3 = p101 + (p111 - p101)*w[6][j] #   5 6
                    edgey4 = p001 + (p011 - p001)*w[7][j] #   4 7

                    edgez1 = p000 + (p001 - p000)*w[8][k] #    0 4
                    edgez2 = p100 + (p101 - p100)*w[9][k] #    1 5
                    edgez3 = p110 + (p111 - p110)*w[10][k]#   2 6
                    edgez4 = p010 + (p011 - p010)*w[11][k]#   3 7


                    # calculate the importance factors for all edges

                    # x-direction
                    impx1 = (
                        (1.0 - w[0][i])*(1.0 - w[4][j])*(1.0 - w[8][k])
                      + w[0][i]*(1.0 - w[5][j])*(1.0 - w[9][k])
                    )

                    impx2 = (
                        (1.0 - w[1][i])*w[4][j]*(1.0 - w[11][k])
                      + w[1][i]*w[5][j]*(1.0 - w[10][k])
                    )

                    impx3 = (
                         (1.0 - w[2][i])*w[7][j]*w[11][k]
                       + w[2][i]*w[6][j]*w[10][k]
                    )

                    impx4 = (
                        (1.0 - w[3][i])*(1.0 - w[7][j])*w[8][k]
                      + w[3][i]*(1.0 - w[6][j])*w[9][k]
                    )

                    magImpx = impx1 + impx2 + impx3 + impx4
                    impx1 /= magImpx
                    impx2 /= magImpx
                    impx3 /= magImpx
                    impx4 /= magImpx


                    # y-direction
                    impy1 = (
                        (1.0 - w[4][j])*(1.0 - w[0][i])*(1.0 - w[8][k])
                      + w[4][j]*(1.0 - w[1][i])*(1.0 - w[11][k])
                    )

                    impy2 = (
                        (1.0 - w[5][j])*w[0][i]*(1.0 - w[9][k])
                      + w[5][j]*w[1][i]*(1.0 - w[10][k])
                    )

                    impy3 = (
                        (1.0 - w[6][j])*w[3][i]*w[9][k]
                      + w[6][j]*w[2][i]*w[10][k]
                    )

                    impy4 = (
                        (1.0 - w[7][j])*(1.0 - w[3][i])*w[8][k]
                      + w[7][j]*(1.0 - w[2][i])*w[11][k]
                    )

                    magImpy = impy1 + impy2 + impy3 + impy4
                    impy1 /= magImpy
                    impy2 /= magImpy
                    impy3 /= magImpy
                    impy4 /= magImpy


                    # z-direction
                    impz1 = (
                        (1.0 - w[8][k])*(1.0 - w[0][i])*(1.0 - w[4][j])
                      + w[8][k]*(1.0 - w[3][i])*(1.0 - w[7][j])
                    )

                    impz2 = (
                        (1.0 - w[9][k])*w[0][i]*(1.0 - w[5][j])
                      + w[9][k]*w[3][i]*(1.0 - w[6][j])
                    )

                    impz3 = (
                        (1.0 - w[10][k])*w[1][i]*w[5][j]
                      + w[10][k]*w[2][i]*w[6][j]
                    )

                    impz4 = (
                        (1.0 - w[11][k])*(1.0 - w[1][i])*w[4][j]
                      + w[11][k]*(1.0 - w[2][i])*w[7][j]
                    )

                    magImpz = impz1 + impz2 + impz3 + impz4
                    impz1 /= magImpz
                    impz2 /= magImpz
                    impz3 /= magImpz
                    impz4 /= magImpz


                    # calculate the correction vectors
                    corx1 = impx1*(p[0][i] - edgex1)
                    corx2 = impx2*(p[1][i] - edgex2)
                    corx3 = impx3*(p[2][i] - edgex3)
                    corx4 = impx4*(p[3][i] - edgex4)

                    cory1 = impy1*(p[4][j] - edgey1)
                    cory2 = impy2*(p[5][j] - edgey2)
                    cory3 = impy3*(p[6][j] - edgey3)
                    cory4 = impy4*(p[7][j] - edgey4)

                    corz1 = impz1*(p[8][k] - edgez1)
                    corz2 = impz2*(p[9][k] - edgez2)
                    corz3 = impz3*(p[10][k] - edgez3)
                    corz4 = impz4*(p[11][k] - edgez4)


                    # multiply by the importance factor

                    # x-direction
                    edgex1 *= impx1
                    edgex2 *= impx2
                    edgex3 *= impx3
                    edgex4 *= impx4

                    # y-direction
                    edgey1 *= impy1
                    edgey2 *= impy2
                    edgey3 *= impy3
                    edgey4 *= impy4

                    # z-direction
                    edgez1 *= impz1
                    edgez2 *= impz2
                    edgez3 *= impz3
                    edgez4 *= impz4


                    # add the contributions
                    obj.data.vertices[vertexNo].co = (
                        edgex1 + edgex2 + edgex3 + edgex4
                      + edgey1 + edgey2 + edgey3 + edgey4
                      + edgez1 + edgez2 + edgez3 + edgez4
                    ) / 3.0

                    obj.data.vertices[vertexNo].co += (
                        corx1 + corx2 + corx3 + corx4
                      + cory1 + cory2 + cory3 + cory4
                      + corz1 + corz2 + corz3 + corz4
                    )
                    vertexNo +=1

    return vertexNo
