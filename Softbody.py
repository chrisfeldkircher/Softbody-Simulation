import pygame
from pygame.locals import *
import random
import numpy as np
from itertools import combinations
import operator

MassPnts = 8
Radius = 40

class Point(object):
    def __init__(self):
        self.pos = np.array((0.0, 0.0)) #position
        self.v = np.array((0.0, 0.0)) # velocity
        self.f = np.array((0.0, 0.0)) # Force
        self.acc = np.array((0.0, 0.0)) #acceleration
     
class Spring(object):
    def __init__(self):
        self.point = np.array((0, 0)) # point indexes
        self.length = np.array((0.0, 0.0)) # rest length
        self.n = np.array((0.0, 0.0)) # normal vector

class obj(object):
    def __init__(self,pnts):
        self.pnts = pnts
        
    # returns a list of all point tupels
    def pointList(self):
        pointList = []
        for point in self.pnts:
            pointList.append((point.pos[0], point.pos[1]))
        return pointList
    
    def renderPoly(self,surfaceObj,color):
        pygame.draw.polygon(surfaceObj, color, self.pointList())

        
    """
    def calcCentroid(self):
        signedArea = 0.
        centroidX,centroidY = 0.,0.
        for idx, pt in enumerate(self.pnts):
            nxt = self.pnts[(idx+1)%len(self.pnts)]
            a = pt.x*nxt.y-pt.y*nxt.x
            signedArea += a
            centroidX += (pt.x+nxt.x)*a
            centroidY += (pt.y+nxt.y)*a
        signedArea *= .5
        return (int(centroidX/(6.*signedArea)),int(centroidY/(6.*signedArea)))
    """
    
class elastic_body(obj):
    def __init__(self,myPoints,mass,final_pressure,ks,kd):
        super(elastic_body,self).__init__(myPoints)
        self.Springs = []         # holds all springs
        self.final_pressure = final_pressure     # max pressure
        self.ks,self.kd = ks,kd     # elasticity(ks) and damping(kd)
        self.mass = mass            # mass of a mass point
        self.width, self.height = 700, 500
        #self.pressure = 10
        self.pressure = 4
        self.dt = 0.08
        self.num_masspnts = MassPnts 
        self.drag = 0.01

        for i in range(0, self.num_masspnts ):
            for j in range(0,i):
                if j != i:
                    self.addSpring(j,i)
    
    def renderPoly(self,surfaceObj):
        super(elastic_body,self).renderPoly(surfaceObj, (109, 22, 138))

    def addSpring(self,i,j):
        # adds a spring between masspoint i and j
        s = Spring()
        s.i, s.j = i,j
        pA,pB = self.pnts[i],self.pnts[j]
        s.length = np.sqrt((pA.pos[0] - pB.pos[0])**2 + (pA.pos[1] - pB.pos[1])**2)
        self.Springs.append(s)
     
    def showSprings(self,surfaceObj):
        for spring in self.Springs:
            p1,p2 = self.pnts[spring.i], self.pnts[spring.j]
            pygame.draw.line(surfaceObj, (21, 214, 140), (p1.pos[0], p1.pos[1]), (p2.pos[0], p2.pos[1]), 2)

    def spring_force(self):
        for s in self.Springs:
            start, end = self.pnts[s.i], self.pnts[s.j]
            x1, y1, x2, y2 = start.pos[0], start.pos[1], end.pos[0], end.pos[1]
        
            dist = np.sqrt((x1-x2)**2 + (y1-y2)**2)
            
            if dist != 0: #start != end
                v_se_x = start.v[0] - end.v[0] # velocity start
                v_se_y = start.v[1] - end.v[1] # velocity end
            
            f = (dist - s.length) * self.ks + ((v_se_x * (x1-x2) + v_se_y * (y1-y2)) * self.kd) / dist #force
            n_x = ((x1-x2) / dist)
            n_y = ((y1-y2) / dist)
            Fx = (n_x * f) # forcevector = n_x*f
            Fy = (n_y * f) # forcevector = n_y*f
            
            start.f[0] -= Fx # set forces for starting point
            start.f[1] -= Fy
            end.f[0] += Fx # set forces for starting point
            end.f[1] += Fy
        
            s.n[0] = (y1-y2) / dist  # calc norm spring_forces
            s.n[1] = (x2-x1) / dist
    
    def pressure_force(self):
        for s in self.Springs: 
            start, end = self.pnts[s.i],self.pnts[s.j]
            x1, y1, x2, y2 = start.pos[0], start.pos[1], end.pos[0], end.pos[1]
            
            dist = np.sqrt((x1-x2)**2 + (y1-y2)**2)
            
            pressure_volume = dist * self.pressure * (1.0/self.calcVolume())
                
            start.f[0] += s.n[0] * pressure_volume
            start.f[1] += s.n[1] * pressure_volume
            end.f[0] += s.n[0] * pressure_volume
            end.f[1] += s.n[1] * pressure_volume

    def update(self, obj, dt):
        new_pos = np.array((obj.pos[0] + obj.v[0]*dt + obj.acc[0]*(dt*dt*0.5), obj.pos[1] + obj.v[1]*dt + obj.acc[1]*(dt*dt*0.5)))
        new_acc = np.array((self.apply_forces(obj)))
        new_v = np.array((obj.v[0] + (obj.acc[0]+new_acc[0])*(dt*0.5), obj.v[1] + (obj.acc[1]+new_acc[1])*(dt*0.5)))
        obj.pos = new_pos;
        obj.v = new_v;
        obj.acc = new_acc;

        if obj.pos[0] > self.width:
            obj.pos[0] = self.width
            obj.v[0] = -obj.v[0]
            
        if obj.pos[1] > self.height:
            obj.pos[1] = self.height
            obj.v[1] = -obj.v[1]
   
    def apply_forces(self, obj):
        grav_acc = [0.0, 9.81 * (abs(self.pressure - self.final_pressure))]
        drag_force = [0.5 * self.drag * (obj.v[0] * abs(obj.v[0])), 0.5 * self.drag * (obj.v[1] * abs(obj.v[1]))]  #D = 0.5 * (rho * C * Area * vel^2)
        drag_acc = [drag_force[0] / self.mass, drag_force[1] / self.mass] #a = F/m
        return (-drag_acc[0]),(grav_acc[1] - drag_acc[1])
    
    def gravity(self):
        for p in self.pnts:
            p.f[0] = 0
            p.f[1] = self.mass * 9.81 * (abs(self.pressure - self.final_pressure))

    def integration(self): #euler-integration - heun methode
        for p in self.pnts:
            drag_force = [0.5 * self.drag * (p.v[0] * abs(p.v[0])), 0.5 * self.drag * (p.v[1] * abs(p.v[1]))]  #D = 0.5 * (rho * C * Area * vel^2)
            drag_acc = [drag_force[0] / self.mass, drag_force[1] / self.mass] #a = F/m

            p.v[0] += ((p.f[0] / self.mass) - drag_acc[0]) * self.dt
            p.pos[0] += p.v[0] * self.dt
            
            #elastic collision
            if p.pos[0] > self.width:
                p.pos[0] = self.width
                p.v[0] = -p.v[0]
            
            if p.pos[0] < 0:
                p.pos[0] = 0
                p.v[0] = -p.v[0]
                    
            p.v[1] += ((p.f[1] / self.mass) - drag_acc[1]) * self.dt
            p.pos[1] += p.v[1] * self.dt
            
            #elastic collision
            if p.pos[1] > self.height:
                p.pos[1] = self.height
                p.v[1] = -p.v[1]
            
            if p.pos[1] < 0:
                p.pos[1] = 0
                p.v[1] = -p.v[1]

    def updatePhysics(self):
        self.gravity()
        self.spring_force()
        self.pressure_force()
        self.integration()

        """
        for p in self.pnts:
            self.update(p, self.dt)
        """
        
    def calcVolume(self):
        self.volume = 0
        # calculate volume --> divergence 
        for s in self.Springs:
            start, end = self.pnts[s.i], self.pnts[s.j]
            x1, y1, x2, y2 = start.pos[0], start.pos[1], end.pos[0], end.pos[1]
        
            dist = np.sqrt((x1-x2)**2 + (y1-y2)**2)
        
            self.volume += .5 * abs(x1-x2) * abs(s.n[0]) * dist
        return self.volume
        
class Softbody(elastic_body):
    def __init__(self,x,y,nump):
        self.x, self.y = x, y
        self.radius = Radius
        # init obj points
        points = []
        for i in range(0,nump):
            p = Point()
            p.pos[0] = self.radius * np.sin(i * (2.0 * np.pi) / nump) + x
            p.pos[1] = self.radius * np.cos(i * (2.0 * np.pi) / nump) + y
            points.append(p)
        #super(Softbody,self).__init__(points,1,3,4,.5) #points, mass, final_pressure, stiffness, damping
        super(Softbody,self).__init__(points,1,6,8,.1) #points, mass, final_pressure, stiffness, damping

class Simulation:
    def __init__(self):
        self.screen = None
        self.dim = self.width, self.height = 700, 500
        self.num_masspnts = MassPnts
        self.obj = []
        self.jump = -180
        self.mx, self.my = 0, 0
        self.apply = False
        self.show_springs = False
            
    def init(self):
        pygame.init()
        pygame.display.set_caption("Softbody Simulation")
        self.screen = pygame.display.set_mode(self.dim)
        self.isRunning = True
        self.screen.fill((255,255,255))
        self.clock = pygame.time.Clock()
        obj1 = Softbody(self.width/2, self.height-Radius, self.num_masspnts)
        self.obj.append(obj1)

    def norm_to_edges(self, polygon):
        unit_norm = []
        for i in range(len(polygon.pnts)):
            next = i + 1
            next %= len(polygon.pnts)

            temp = [polygon.pnts[i].pos[0] - polygon.pnts[next].pos[0], polygon.pnts[i].pos[1] - polygon.pnts[next].pos[1]]
            norm = np.array((temp[1], -temp[0]))
            unit_n = np.array((norm[0]/np.linalg.norm(norm), norm[1]/np.linalg.norm(norm)))
            unit_norm.append(tuple(unit_n))
        return unit_norm

    def coillision(self, polygon1, polygon2):
        #seperate axis theorem --> works only with convex shapes
        #1 get the normals to each side
        #2 project the max and min intervalls on the normals
        #3 see if the intervalls overlap --> all have to overlap for a collision
        #4 calculate the minimal push vector to seperate the polygons

        axis1 = []
        axis2 = []
        push = []
        min1, max1 = float('+inf'), float('-inf')
        min2, max2 = float('+inf'), float('-inf')
        normals_to_test = [norm for norm in self.norm_to_edges(polygon1)]
        pnts1 = []
        pnts2 = []

        for i in range(len(polygon1.pnts)):
            pnts1.append(polygon1.pnts[i].pos)
            pnts2.append(polygon2.pnts[i].pos)

        for orth in normals_to_test:
            sep, pushvec = self.proj_on_axis(orth, pnts1, pnts2)

            if not sep:
                push.append(pushvec) 
            else:
                return False

        mpv =  min(push, key=(lambda v: np.dot(v, v)))
        d = self.centers_displacement(pnts1, pnts2)  
        if np.dot(d, mpv) > 0: # if it's the same direction, then invert
            mpv[0] *= -1
            mpv[1] *= -1

        for p in polygon1.pnts:
            p.v[0] += mpv[0]
            p.v[1] += mpv[1]

    def proj_on_axis(self, norm, pnts1, pnts2):
        min1, max1 = float('+inf'), float('-inf')
        min2, max2 = float('+inf'), float('-inf')

        for vert in pnts1:
            proj = np.dot(vert, norm)

            min1 = min(min1, proj)
            max1 = max(max1, proj)

        for vert in pnts2:
            proj = np.dot(vert, norm)

            min2 = min(min2, proj)
            max2 = max(max2, proj)

        if max1 >= min2 and max2 >= min1:
            d = min(max2 - min1, max1 - min2)
            d_over_norm_squared = d/np.dot(norm, norm) + 1e-10
            pushvec = [d_over_norm_squared*norm[0], d_over_norm_squared*norm[1]]
            return False, pushvec
        else:
            return True, None 
    
    def centers_displacement(self, p1, p2):
        # geometric center
        c1 = np.mean(np.array(p1), axis=0)
        c2 = np.mean(np.array(p2), axis=0)
        return c2 - c1

    def apply_force(self):
        """
        if self.apply:
            self.mx, self.my = pygame.mouse.get_pos()
            #force to the first three masspoints in mouse pointer direction
            for obj in self.obj:
                pnt1 = obj.pnts[1]
                pnt2 = obj.pnts[0]
                pnt3 = obj.pnts[2]
                x = self.mx - pnt1.pos[0]
                y = self.my - pnt1.pos[1]
                pnt1.v[0], pnt1.v[1] =  x, y
                pnt2.v[1], pnt3.v[1], pnt2.v[0], pnt3.v[0] = y/2, y/2, x/2, x/2
        """
        if self.apply:
            self.mx, self.my = pygame.mouse.get_pos()
            #force to the first two masspoints in mouse pointer direction
            for obj in self.obj:
                pnt1 = obj.pnts[1]
                pnt2 = obj.pnts[0]
                x = self.mx - pnt1.pos[0]
                y = self.my - pnt1.pos[1]
                pnt1.v[0], pnt1.v[1] =  x, y
                pnt2.v[1], pnt2.v[0], = y/2, x/2

    def event(self, event):
        if event.type == pygame.QUIT:
            self.isRunning = False
        elif event.type == pygame.MOUSEBUTTONDOWN:
            if event.button == 1:
                self.apply = True
        elif event.type == pygame.MOUSEBUTTONUP:
            if event.button == 1:
                self.apply = False

        elif event.type == pygame.KEYDOWN:
            if event.key == pygame.K_SPACE:
                for obj in self.obj:
                    for p in obj.pnts:
                        p.v[1] += self.jump
            if event.key == pygame.K_c:
                self.mx, self.my = pygame.mouse.get_pos()
                self.obj.append(Softbody(self.mx, self.my, self.num_masspnts))
            if event.key == K_a:
                self.show_springs = not self.show_springs
            if event.key == K_LEFT:
               for obj in self.obj:
                    for p in obj.pnts:
                        p.v[0] -= 10
            if event.key == K_RIGHT:
                   for obj in self.obj:
                    for p in obj.pnts:
                        p.v[0] += 10
            if event.key == K_DOWN:
                   for obj in self.obj:
                    for p in obj.pnts:
                        p.v[1] += 10

    def loop(self):
        pairs = combinations(range(len(self.obj)), 2)
        for i,j in pairs:
            self.coillision(self.obj[i], self.obj[j])

        for obj in self.obj:
            obj.updatePhysics()
            pass
        self.apply_force()

    def draw(self):
        self.screen.fill((255,255,255))

        for obj in self.obj:
            obj.renderPoly(self.screen)
            if self.show_springs:
                obj.showSprings(self.screen)
            #pygame.draw.circle(self.screen, (0,100,0), obj.centers_displacement(), 5, 3)
        pygame.display.flip()

    def execute(self):
        if self.init() == False:
            self.isRunning = False

        while self.isRunning:
            for event in pygame.event.get():
                self.event(event)
            
            #self.clock.tick(40)
            self.loop()
            self.draw()

        pygame.quit()

if __name__ == "__main__":
    p = Simulation()
    p.execute()
