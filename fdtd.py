from __future__ import division

import sys
import signal
import math

import numpy
import pygame

pygame.init()
signal.signal(signal.SIGQUIT, signal.SIG_DFL)

SIZE = 320, 240

Ex = numpy.zeros(SIZE)
Ey = numpy.zeros(SIZE)
if 0:
    for x in xrange(SIZE[0]):
        for y in xrange(SIZE[1]):
            dx = x - (SIZE[0]-1)/2
            dy = y - (SIZE[1]-1)/2
            dist = math.sqrt(dx**2 + dy**2)
            Ex[x, y] = 1e3*dx/dist**2
            Ey[x, y] = 1e3*dy/dist**2
Hz = numpy.zeros(SIZE)
notmetal = numpy.ones(SIZE, dtype=numpy.uint8)
er = numpy.ones(SIZE)
sigma = numpy.zeros(SIZE)
batteries = []

def reset():
    er.fill(1)
    notmetal.fill(1)
    sigma.fill(0)
    ML_THICKNESS = 25
    for i in xrange(ML_THICKNESS):
        k = i/ML_THICKNESS
        s = k**3 * .2
        sigma[ML_THICKNESS-1-i, :] = s
        sigma[-(ML_THICKNESS-1-i)-1, :] = s
        sigma[:, ML_THICKNESS-1-i] = s
        sigma[:, -(ML_THICKNESS-1-i)-1] = s
    batteries[:] = []
reset()

d = pygame.display.set_mode(SIZE)

def draw(pos):
    for dx in xrange(-5, 5+1):
        for dy in xrange(-5, 5+1):
            if math.hypot(dx, dy) <= 5:
                Hz[pos[0]+dx, pos[1]+dy] += 20
def draw_metal(pos, erase):
    for dx in xrange(-5, 5+1):
        for dy in xrange(-5, 5+1):
            notmetal[pos[0]+dx, pos[1]+dy] = erase
            if not erase:
                Hz[pos[0]+dx, pos[1]+dy] = 0
def draw_dielectric(pos):
    for dx in xrange(-5, 5+1):
        for dy in xrange(-5, 5+1):
            er[pos[0]+dx, pos[1]+dy] = 4
            sigma[pos[0]+dx, pos[1]+dy] = .02
def draw_battery((pos, inv, mag)):
    if mag:
        for dx in xrange(-5, 5+1):
            for dy in xrange(-5, 5+1):
                Hz[pos[0]+dx, pos[1]+dy] = 300*math.sin(2*math.pi*20*t) * (-1 if inv else 1)
        return
    for dx in xrange(-5, 5+1):
        for dy in xrange(-5, 5+1):
            Ex[pos[0]+dx, pos[1]+dy] = 0
            Ey[pos[0]+dx, pos[1]+dy] = 1000 * (-1 if inv else 1) #* math.sin(2*math.pi*20*t)

pixel_size = 1e-3


ts = 0
t = 0
keys = set()
pause = False
mouse_down_pos = None
while True:
    for event in pygame.event.get():
        #print event
        
        if event.type == pygame.KEYDOWN: keys.add(event.key)
        if event.type == pygame.KEYUP: keys.discard(event.key)
        
        if event.type == pygame.QUIT:
            sys.exit()
        elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
            mouse_down_pos = event.pos
            #draw(event.pos)
        elif event.type == pygame.MOUSEBUTTONUP and event.button == 1:
            mouse_down_pos = None
            #draw(event.pos)
        elif event.type == pygame.MOUSEMOTION and event.buttons[0]:
            #draw(event.pos)
            mouse_down_pos = event.pos
        elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 3:
            draw_metal(event.pos, pygame.K_LSHIFT in keys)
        elif event.type == pygame.MOUSEMOTION and event.buttons[2]:
            draw_metal(event.pos, pygame.K_LSHIFT in keys)
        elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 2:
            draw_dielectric(event.pos)
        elif event.type == pygame.MOUSEMOTION and event.buttons[1]:
            draw_dielectric(event.pos)
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_c:
            print event
            Ex.fill(0)
            Ey.fill(0)
            Hz.fill(0)
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_r:
            reset()
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_b:
            batteries.append((pygame.mouse.get_pos(), event.mod & pygame.KMOD_SHIFT, event.mod & pygame.KMOD_CTRL))
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_p:
            pause = not pause
    
    dt = 3e-4
    
    if not pause:
        Jx = numpy.zeros((SIZE[0], SIZE[1]-1))
        Jy = numpy.zeros((SIZE[0]-1, SIZE[1]))
        if mouse_down_pos is not None:
            for dx in xrange(-5, 5+1):
                for dy in xrange(-5, 5+1):
                    Jx[mouse_down_pos[0]+dx, mouse_down_pos[1]+dy] = 1e6
                    Jy[mouse_down_pos[0]+dx, mouse_down_pos[1]+dy] = 0
        for i in xrange(5):
            if True: #if ts % 2 == 0:
                curlHx = (-Hz[:, :-1] + Hz[:, 1:])/pixel_size # defined for [:, :-1]
                curlHy = (Hz[:-1, :] - Hz[1:, :])/pixel_size # defined for [:-1, :]
                #Ex[:, :-1] += dt * curlHx / er[:, :-1]
                #Ey[:-1, :] += dt * curlHy / er[:-1, :]
                #Ex[:, :-1] += dt * curlHx
                #Ey[:-1, :] += dt * curlHy
                # dE/dt = (curl H - J - sigma E)/er
                curlHx -= Jx
                curlHy -= Jy
                # equations continue without J
                # newE = oldE + dt * ((curl H - sigma E)/er)
                # newE = oldE + dt * ((curl H - sigma (newE+oldE)/2)/er)
                # (1 + sigma/2 er) newE = oldE + dt * curl H / er - sigma/2er oldE
                # newE = (dt * curl H / er + (1 - sigma/2er) oldE)/(1 + sigma/2 er)
                Ex[:, :-1] = ((1-sigma[:, :-1]/2/er[:, :-1]) * Ex[:, :-1] + dt*curlHx/er[:, :-1])/(1 + sigma[:, :-1]/2/er[:, :-1])
                Ey[:-1, :] = ((1-sigma[:-1, :]/2/er[:-1, :]) * Ey[:-1, :] + dt*curlHy/er[:-1, :])/(1 + sigma[:-1, :]/2/er[:-1, :])
                map(draw_battery, batteries)
                Ex *= notmetal
                Ey *= notmetal
                t += dt/2
            if True:#else:
                curlE = (Ey[1:, 1:] - Ex[1:, 1:] - Ey[:-1, 1:] + Ex[1:, :-1])/pixel_size # defined for [1:, 1:]
                # dH/dt = -curl E / mu_r
                Hz[1:, 1:] -= dt * curlE #/ mu_r
                t += dt/2
            ts += 2#1
    
    
    if pygame.mouse.get_pressed()[0]:
        draw(pygame.mouse.get_pos())
    
    print numpy.min(Ex), numpy.max(Ex)
    
    pygame.surfarray.pixels3d(d)[:, :, 0] = numpy.round(128+Ex)
    pygame.surfarray.pixels3d(d)[:, :, 1] = numpy.round(128+Ey)
    pygame.surfarray.pixels3d(d)[:, :, 2] = numpy.round(128+Hz)
    pygame.surfarray.pixels3d(d)[:, :, 0] *= notmetal
    pygame.surfarray.pixels3d(d)[:, :, 1] *= notmetal
    pygame.surfarray.pixels3d(d)[:, :, 0] += (10*er).astype(numpy.uint8)
    pygame.surfarray.pixels3d(d)[:, :, 1] += (1000*sigma).astype(numpy.uint8)
    pygame.display.update()
