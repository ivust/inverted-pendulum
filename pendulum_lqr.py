from math import *
from numpy import *
import control.matlab
import pygame

dt = 0.01
g = 9.81
l = 1.0
m = 1.0

global q1, q2, q3, q4, r
q1 = 0.001
q2 = 0.001
q3 = 100000
q4 = 2000
r = 0.005



clock = pygame.time.Clock()
pygame.init()
size = (width, height) = (1200,600)
screen = pygame.display.set_mode(size)

class Pendulum:
	def __init__(self, x0, y0, angle0):
		self.angle0 = angle0
		self.angle = angle0
		self.velocity = 0
		self.x0 = x0
		self.y0 = y0
		self.x0_vel = 0
		self.x = x0 + 250.0 * sin(angle0)
		self.y = y0 + 250.0 * cos(angle0)
	def move(self, control):
		self.angle = atan2(self.x - self.x0, self.y - self.y0)
		d_velocity = -g * sin(self.angle) * dt / l
		self.velocity = self.velocity + d_velocity
		d_angle = dt * self.velocity
		self.angle = self.angle + d_angle
		self.x = self.x0 + 250.0 * sin(self.angle)
		self.y = self.y0 + 250.0 * cos(self.angle)

		d_x0_vel = dt * control
		self.x0_vel = self.x0_vel + d_x0_vel
		dx0 = dt * self.x0_vel
		self.x0 = self.x0 + dx0

		if self.x0 > 1200 or self.x0 < 0:
			self.x0 = self.x0 - dx0

	def draw(self):
		pygame.draw.circle(screen, (0,0,0), [int(self.x0), int(self.y0)], 5)
		pygame.draw.line(screen, (0,0,0), [self.x0, self.y0], [self.x, self.y], 2)
		pygame.draw.circle(screen, (255,0,0), [int(self.x), int(self.y)], 10)
		pygame.draw.line(screen, (0,0,0), [0, self.y0], [1200, self.y0], 3)

class LQR:
	def __init__(self, q1, q2, q3, q4, r, pendulum):
		self.q1 = q1
		self.q2 = q2
		self.q3 = q3
		self.q4 = q4
		self.p = pendulum
		self.A = matrix([[0,1,0,0],[0,0,-g,0],[0,0,0,1],[0,0,2*g,0]])
		self.B = matrix([[0],[1],[0],[-1]])
		self.Q = diag([q1,q2,q3,q4])
		self.R = r
		self.K = control.matlab.lqr(self.A, self.B, self.Q, self.R)[0]
		print self.K
	def update(self):
		self.q1 = q1
		self.q2 = q2
		self.q3 = q3
		self.q4 = q4
		self.Q = diag([q1,q2,q3,q4])
		self.K = control.matlab.lqr(self.A, self.B, self.Q, self.R)[0]

		
	def output(self):
		# X = matrix([[self.p.x0 - 510], [self.p.x0_vel],
		# 		[copysign(1, self.p.angle) * (pi - abs(self.p.angle))], [self.p.velocity]])
		X = matrix([[-(self.p.x0 - 600) / 10], [self.p.x0_vel / 10],
				[copysign(1, self.p.angle) * (-pi + abs(self.p.angle))], [self.p.velocity]])
		U = self.K * X
		# print U.flat[0], self.p.velocity, self.p.x0_vel, self.p.x0 - 500
		# print self.K
		return U.flat[0]
		# return 0

def lqr_tune(lqr, pend):
	global q1, q2, q3, q4
	n_params = 4
	init_params = [lqr.q1, lqr.q2, lqr.q3, lqr.q4]
	dparams    = [0.0001, 0.0001, 1000, 10]
	params     = [0.0 for i in range(n_params)]

	for i in range(n_params):
		params[i] = init_params[i]

	best_error = 0.0;
	for t in range(100):
		screen.fill((255,255,255))
		lqr.update()
		pend.move(lqr.output())
		best_error = best_error + abs(pend.angle)
		pend.draw()
		clock.tick(60)
		pygame.display.flip()
	best_error = best_error / 100
	best_error = 3 * (pi - best_error)
	print best_error

	n = 1
	while sum(dparams) > 500 and n <= 15:
		# print dparams
		for i in range(len(params)):
			params[i] += dparams[i]
			err = 0

			[q1, q2, q3, q4] = params
			for t in range(100):
				screen.fill((255,255,255))
				lqr.update()
				pend.move(lqr.output())
				err = err + abs(pend.angle)
				pend.draw()
				clock.tick(60)
				pygame.display.flip()

			err = err / 100
			err = pi - err
			err = (1 + 1/n)*err


			if err < best_error:
				best_error = err
				dparams[i] *= 1.1
			else:
				params[i] -= 2.0 * dparams[i]            
				err = 0
				
				[q1, q2, q3, q4] = params
				for t in range(100):
					screen.fill((255,255,255))
					lqr.update()
					pend.move(lqr.output())
					err = err + abs(pend.angle)
					pend.draw()
					clock.tick(60)
					pygame.display.flip()

				err = err / 100
				err = pi - err
				err = (1 + 1/n)*err

				if err < best_error:
					best_error = err
					dparams[i] *= 1.1
				else:
					params[i] += dparams[i]
					dparams[i] *= 0.5
		n += 1
		print '#', n, params, ' -> ', best_error
		print dparams
	print ' '

	return params


def draw_text():
	myfont = pygame.font.SysFont("monospace", 15)
	label1 = myfont.render("Cart position parameter: %.5f" % q1, 1, (255,0,0))
	screen.blit(label1, (50, 500))
	label2 = myfont.render("Cart velocity parameter: %.5f" % q2, 1, (255,0,0))
	screen.blit(label2, (50, 520))
	label3 = myfont.render("Angular positon parameter: %.1f" % q3, 1, (255,0,0))
	screen.blit(label3, (50, 540))
	label3 = myfont.render("Angular velocity parameter: %.1f" % q4, 1, (255,0,0))
	screen.blit(label3, (50, 560))

def make_buttons(q1, q2, q3, q4, lqr, pend):
	pygame.draw.rect(screen, (0, 0, 255), [420, 500, 90, 15])
	pygame.draw.rect(screen, (0, 0, 255), [520, 500, 90, 15])
	pygame.draw.rect(screen, (0, 0, 255), [420, 520, 90, 15])
	pygame.draw.rect(screen, (0, 0, 255), [520, 520, 90, 15])
	pygame.draw.rect(screen, (0, 0, 255), [420, 540, 90, 15])
	pygame.draw.rect(screen, (0, 0, 255), [520, 540, 90, 15])
	pygame.draw.rect(screen, (0, 0, 255), [420, 560, 90, 15])
	pygame.draw.rect(screen, (0, 0, 255), [520, 560, 90, 15])
	pygame.draw.rect(screen, (0, 0, 255), [700, 500, 120, 55])

	myfont = pygame.font.SysFont("monospace", 15)
	label1 = myfont.render("Increase", 1, (255,255,255))
	label2 = myfont.render("Decrease", 1, (255,255,255))
	label3 = myfont.render("Tune", 1, (255,255,255))
	label4 = myfont.render("automatically", 1, (255,255,255))
	screen.blit(label1, (420, 500))
	screen.blit(label2, (520, 500))
	screen.blit(label1, (420, 520))
	screen.blit(label2, (520, 520))
	screen.blit(label1, (420, 540))
	screen.blit(label2, (520, 540))
	screen.blit(label1, (420, 560))
	screen.blit(label2, (520, 560))
	screen.blit(label3, (740, 500))
	screen.blit(label4, (702, 530))

	if(pygame.mouse.get_pressed()[0]):
		pos = (pos1, pos2) = pygame.mouse.get_pos()
		if (420 <= pos1 <= 510 and 500 <= pos2 <= 515):
			q1 = q1 + 0.001
		elif (520 <= pos1 <= 590 and 500 <= pos2 <= 515):
			q1 = q1 - 0.001
			if (q1 < 0.001): q1 = q1 + 0.001
		elif (420 <= pos1 <= 510 and 520 <= pos2 <= 535):
			q2 = q2 + 0.001
		elif (520 <= pos1 <= 590 and 520 <= pos2 <= 535):
			q2 = q2 - 0.001
		elif (420 <= pos1 <= 510 and 540 <= pos2 <= 555):
			q3 = q3 + 1000
		elif (520 <= pos1 <= 590 and 540 <= pos2 <= 555):
			q3 = q3 - 1000
		elif (420 <= pos1 <= 510 and 560 <= pos2 <= 575):
			q4 = q4 + 10
		elif (520 <= pos1 <= 590 and 560 <= pos2 <= 575):
			q4 = q4 - 10
		elif (700 <= pos1 <= 820 and 500 <= pos2 <= 555):
			(q1, q2, q3, q4) = lqr_tune(lqr, pend)

	return (q1, q2, q3, q4)



p = Pendulum(500, 300, pi - pi/8)
lqr = LQR(q1, q2, q3, q4, r, p)
while 1:
	screen.fill((255,255,255))
	for event in pygame.event.get():
		if event.type == pygame.QUIT:
			pygame.quit()
	lqr.update()
	p.move(lqr.output())
	# p.move(10)
	p.draw()
	draw_text()
	(q1, q2, q3, q4) = make_buttons(q1, q2, q3, q4, lqr, p)


	clock.tick(60)
	pygame.display.flip()