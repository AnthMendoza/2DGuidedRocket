import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import csv



gavitationalAcceleration = 9.8 #m/s
airDensity = 1.225 #kg/m^3
positionVector = [0,500,0] # m [x,y]
vehicleState = [0,1,0] # 
velocityVector = [0,0,0] #m/s [x,y]
airVelocityVector = [0,0,0] # m/s [x,y]
acclerationVector = [0,0,0] # m/s^2 [x,y]
forceVector = [0,0,0] # newtons [x,y]
targetPosition = [80,0,0] #m [x,y]
timeStep = .005 # s
splineIndex = 0


logsVehicleVector=[]
logsXPosition=[]
logsYPosition=[]
logsLift = []
logsDrag=[]
logsVelocity = []
logsVehicleVector=[]
logsVehicleAirVelocity = []
logsAppliedFinForce = []
logsRotation = []
logsTimeStamp =  []

vehicleMass = 1.5 #kg
vehicleMOI =[.03,.03,0.1] #[pitch, yaw, roll]  
rotatinoalVelocity = 0

vehicleCenterOfPressureConstant = [0,.25,0]#m [x,y,z] from center rear of the vehicle, center axis of the body
rotationalDampingFactor = 10
vehicleCenterOfGravity = [0,.25,0]#m [x,y,z] from center rear of the vehicle, center axis of the body
bezierSplinePoints = []

csvFileName = f'sim_S{positionVector}T{targetPosition}W{airVelocityVector}V{velocityVector}.csv'
#need to compute new desired minimum trijectory 


def preFlight():
    global bezierSplinePoints
    bezierSplinePoints = bezierSpline(targetPosition,positionVector)


def mainLoop():
    preFlight()
    global iterations
    iterations = 0
    while positionVector[1] > targetPosition[1]:
        iterations = iterations + 1
        sumof = sumOfMomentsAndAccels()
        velocityVectorUpdate = velocityUpdate(sumof[1],timeStep)
        global vehicleState
        vehicleState = rotationalPositionUpdate(rotatinoalVelocityUpdate(rotationalAccelUpdate(sumof[0]),timeStep),timeStep)
        

        logsRotation.append(np.degrees(np.arctan(vehicleState[0]/vehicleState[1])))
        

        for v in range(len(velocityVector)):
            
            velocityVector[v] = velocityVector[v] + velocityVectorUpdate[v]

        

        deltaPosition = positionUpdate(velocityVector , timeStep)
        for i in range(len(positionVector)):
            positionVector[i] = positionVector[i] + deltaPosition[i]
            indexVector = [0,1,0]
            #theta = angleOfTwoVectors(indexVector,vehicleState)
        logsVelocity.append(np.linalg.norm(velocityVector))
        logsVehicleAirVelocity.append(np.add(velocityVector , airVelocityVector))
        logsVehicleVector.append(vehicleState)
        logsXPosition.append(positionVector[0])
        logsYPosition.append(positionVector[1])
        logsTimeStamp.append(iterations * timeStep)
    print(f"Time to hit the ground in seconds {timeStep*iterations} number of iterations {iterations} landed at {logsXPosition[-1]}")
        

        
   

        

        
def relativeAirSpeed():
    relativeAirSpeed = []
    for i , velo in enumerate(velocityVector):
        relativeAirSpeed.append(velo + airVelocityVector[i])
    return relativeAirSpeed
    


def angleOfTwoVectors(vector1,vector2):
    vect1 = np.array(vector1)
    vect2 = np.array(vector2)
    
    if np.linalg.norm(vect1) == 0 or np.linalg.norm(vect2) == 0:
        return 0
    cross = np.cross(vector1,vector2)

    if cross[2] < 0:
        return -np.arccos((np.dot(vect1,vect2))/(np.linalg.norm(vect1)*np.linalg.norm(vect2)))
    return np.arccos((np.dot(vect1,vect2))/(np.linalg.norm(vect1)*np.linalg.norm(vect2)))





def vehicleCenterOfPressure(theta,const):   #Calling this function due to real world Center of pressure will change with alpha angle of attack
    return vehicleCenterOfPressureConstant  #m from rear center of the vehicle



def dragConstants(theta):   
    area = abs(np.sin(theta))*0.0258064+.00811
    return (abs(np.sin(theta ))*.5+.55)*area*.5



def liftConstants(theta):
    area = abs(np.sin(theta*2))*0.0258064#+.00811
    return(abs(np.sin(theta*2))*1+.2)*area



def dragForce(LOCvelocityVector,LOCairVelocityVector,LOCvehicleState): #tested 
    adjustedVelocityVector = np.add(LOCvelocityVector , LOCairVelocityVector)
    theta = angleOfTwoVectors(adjustedVelocityVector,LOCvehicleState)
    absoluteVelocity= np.linalg.norm(adjustedVelocityVector)
    if absoluteVelocity == 0:
        normalized = [0,0,0]
    else:
        normalized = adjustedVelocityVector/absoluteVelocity
    
    normalized = np.array(normalized)

    dragVector = normalized * (1/2 *(airDensity * dragConstants(theta) * absoluteVelocity**2 ))
    dragVector[0] = - dragVector[0]
    dragVector[1] = - dragVector[1]
    logsDrag.append(dragVector)

    return dragVector



def liftForce(LOCvelocityVector,LOCairVelocityVector,LOCvehicleState): #tested 
    adjustedVelocityVector = -np.add(LOCvelocityVector , LOCairVelocityVector)
    adjustedVelocityVector[0] = -adjustedVelocityVector[0]
    #adjustedVelocityVector[1] = adjustedVelocityVector[1]
    theta = angleOfTwoVectors(adjustedVelocityVector,LOCvehicleState)
    absoluteVelocity= np.linalg.norm(adjustedVelocityVector)
    crossProduct = np.cross(adjustedVelocityVector,LOCvehicleState)
    if crossProduct[2] > 0:
        perpendicularVector = np.array([-adjustedVelocityVector[1],-adjustedVelocityVector[0],0])
        normilzedperpendicularVector = perpendicularVector/np.linalg.norm(perpendicularVector)
    elif crossProduct[2]<0:
        perpendicularVector = np.array([adjustedVelocityVector[1],adjustedVelocityVector[0],0])
        normilzedperpendicularVector = perpendicularVector/np.linalg.norm(perpendicularVector)
    else:
        normilzedperpendicularVector = np.array([0,0,0])
    
    liftVector = normilzedperpendicularVector* (1/2 *(airDensity * liftConstants(theta) * absoluteVelocity**2 ))

    logsLift.append(liftVector)
    return liftVector



def rotationalDamping(currentRotationalVelocity,MOI,dampingFactor): #rotaional Damping force will need to be estimated based on simulated/empirical data
    adjustmentForce=-(currentRotationalVelocity)*dampingFactor
    return adjustmentForce*MOI




def finForce(LOCvehicleState): 
    orthogonalVector = [-LOCvehicleState[1],LOCvehicleState[0],0]
    orthogonalVectorLength = np.linalg.norm(orthogonalVector)
    if orthogonalVectorLength == 0:
        normalized = [0,0,0]
    else:
        normalized = orthogonalVector/orthogonalVectorLength
    normalized = np.array(normalized)
    return normalized



def sumOfMomentsAndAccels():
    moment = 0
    result1 = forceApplied(dragForce(velocityVector,airVelocityVector,vehicleState),vehicleState,vehicleCenterOfPressureConstant,vehicleMass,vehicleMOI[0])

    result2 = forceApplied(liftForce(velocityVector,airVelocityVector,vehicleState),vehicleState,vehicleCenterOfPressureConstant,vehicleMass,vehicleMOI[0])

    stanleyOuput = Stanley.compute(velocityVector,airVelocityVector,positionVector,bezierSplinePoints[bezierIndex(positionVector,splineIndex,bezierSplinePoints)])
    PDOutput = controller.compute(stanleyOuput,angleOfTwoVectors([0,1,0],vehicleState),iterations*timeStep)
    PDOutput[0] = -PDOutput[0]
    result3 = forceApplied(PDOutput,vehicleState,vehicleCenterOfGravity,vehicleMass,vehicleMOI[0])
    
    moment = moment  +result1[0] -result2[0]  + rotationalDamping(rotatinoalVelocity,vehicleMOI[0],rotationalDampingFactor) - result3[0]
    
    
    logsAppliedFinForce.append(result3[1])

    
    npresult1 = np.array(result1[1])
    npresult2 = np.array(result2[1])
    npresult3 = np.array(result3[1])
    force =  [0,-gavitationalAcceleration*vehicleMass,0]+ npresult1  +npresult2 #+ npresult3


    accelVector = force/vehicleMass
    momentAndAccel = [moment,accelVector]
    return momentAndAccel




def forceApplied(forceVector,vehicleVector,forceIncident,mass,MOI): #tested, tested in edge case where vector lenght = 0

    npVehicleVector = np.array(vehicleVector)
    lengthOfVector = np.linalg.norm(npVehicleVector)
    npForceVector = np.array(forceVector)
    npForceIncident = np.array(forceIncident)

    vehicleVectorNormilized = npVehicleVector * np.linalg.norm(vehicleVector)
    
    vehicleVectorNormilized = vehicleVectorNormilized * forceIncident[1]
    crossProduct = np.cross(npForceVector,vehicleVectorNormilized)


    for v, vect in enumerate(npVehicleVector):
        if lengthOfVector > 0:
            npVehicleVector[v] = vect/lengthOfVector
        else:
            npVehicleVector[v] = 0
        npVehicleVector * npForceIncident[1]


    if crossProduct[2]> 0:
        return [-(np.linalg.norm(crossProduct)/MOI),forceVector]
    else:
        return [(np.linalg.norm(crossProduct)/MOI),forceVector]


def positionUpdate(velocity,timeStep):
    deltaPosition=[]
    for  velo in velocity:
        deltaPosition.append(velo * timeStep)
    return deltaPosition


#print(forceApplied(dragForce([0,-1,0],airVelocityVector,[1,0,0]),[1,0,0],vehicleCenterOfPressureConstant,vehicleMass,vehicleMOI[0]))
#print(dragForce([0,-1,0],airVelocityVector,[1,0,0]))
#print(forceApplied(liftForce([0,-1,0],airVelocityVector,[1,1,0]),[1,1,0],vehicleCenterOfPressureConstant,vehicleMass,vehicleMOI[0]))


def velocityUpdate(accelVector,timeStep):
    deltaVelo=[]
    for  accel in accelVector:
        deltaVelo.append(accel * timeStep)
    return deltaVelo



def rotationalAccelUpdate(moment):
    return moment/vehicleMOI[0]


def rotatinoalVelocityUpdate(accel,timestep):
    global rotatinoalVelocity
    rotatinoalVelocity = rotatinoalVelocity + accel*timestep
    
    return rotatinoalVelocity



def rotationalPositionUpdate(rotationalVelocity,timestep):
    deltaRadians = rotationalVelocity*timestep
    currentRadians=angleOfTwoVectors([0,1,0],vehicleState)
    currentRadians = currentRadians + deltaRadians
    return [-np.sin(currentRadians),np.cos(currentRadians),0]



def maxFinForce(LOCvelocityVector,LOCairVelocityVector,LOCvehicleState):
    adjustedVelocityVector = np.add(LOCvelocityVector , LOCairVelocityVector)
    adjustedVelocityVector[1] = -adjustedVelocityVector[1]
    theta = angleOfTwoVectors(adjustedVelocityVector,LOCvehicleState)

    vehicleVelocityHeading=abs(np.cos(theta)*np.linalg.norm(adjustedVelocityVector))
    maxForce = .01*(vehicleVelocityHeading**2)
    
    return maxForce





class PDController:
    def __init__(self, Kp, Kd):
        self.Kp = Kp
        self.Kd = Kd
        self.previousError = 0
        self.previousTime = None

    def compute(self, setpoint, measuredValue, currentTime):

        error = setpoint - measuredValue


        if self.previousTime is None:
            dt = 0
        else:
            dt = currentTime - self.previousTime


        if dt > 0:
            derivative = (error - self.previousError) / dt
        else:
            derivative = 0


        force = self.Kp * error + self.Kd * derivative

        self.previousError = error
        self.previousTime = currentTime
        maxForce = maxFinForce(velocityVector,airVelocityVector,vehicleState)
        forceVector  = finForce(vehicleState)

        if abs(force) > maxForce:
            if force < 0:
                
                return -maxForce * forceVector
            else:
                
                return maxForce * forceVector
        
        return force * forceVector




class StanleyController:
    def __init__(self,gainK):
        self.gainK = gainK
        self.lastPointOfIntrest = [0,0,0]
    
    def compute(self, velocityVector,airVelocityVector,currentPosition, currentPointOfIntrest):
        adjustedVelocityVector = np.add(velocityVector , airVelocityVector)
        adjustedVelocity = np.linalg.norm(adjustedVelocityVector)

        deltaX =   currentPointOfIntrest[0]  - self.lastPointOfIntrest[0]

        deltaY =  currentPointOfIntrest[1] -self.lastPointOfIntrest[1]

        if deltaX == 0:
            slope = 0
        else:
            slope = deltaY/deltaX

        crossTrackError = currentPointOfIntrest[0]- currentPosition[0]

        phi = np.arctan(slope)
        

        theta = angleOfTwoVectors([1,0,0],adjustedVelocityVector)
        
        deltaAngle = theta - phi + np.pi/2
        self.lastPointOfIntrest = currentPointOfIntrest
        if adjustedVelocity == 0 :
            output = deltaAngle
        else:
            output =   np.arctan( self.gainK * crossTrackError/adjustedVelocity)
        #print(f" {deltaAngle} {np.arctan( self.gainK * crossTrackError/adjustedVelocity)} {crossTrackError}")
        
        return -output


        




controller = PDController(Kp=40, Kd=.8)

Stanley = StanleyController(gainK=10)





#def controlLoop(forceVector,currentPosition):
#    pGain = 1.5
#    global splineIndex
#    while currentPosition[1] < bezierSplinePoints[splineIndex][1] and splineIndex < len(bezierSplinePoints)-1:
#        
#        splineIndex = splineIndex + 1
#    delta = bezierSplinePoints[splineIndex][0]-currentPosition[0]
#    force  = -delta * pGain
#    maxForce = maxFinForce(velocityVector,airVelocityVector,vehicleState)
#
#    if abs(force) > maxForce:
#        if force < 0:
#            #logsAppliedFinForce.append(-maxForce)
#            return -maxForce * forceVector
#        else:
#            #logsAppliedFinForce.append(maxForce)
#            return maxForce * forceVector
#    #logsAppliedFinForce.append(force)
#    return force * forceVector



def bezierSpline(target,current):
    numOfsteps = 1000
    steps = np.linspace(0,1,numOfsteps)
    xZero = current[0]
    yZero = current[1]
    xOne = current[0]
    yOne = (current[1]-target[1])*1/2
    xTwo = target[0]
    yTwo = (current[1]-target[1])/5
    xThree = target[0]
    yThree = target[1]
    spline = []
    for t in steps:
        x = (1-t)*((1-t)*((1-t)*xZero + t*xOne) + t*((1-t)*xOne+t*xTwo)) + t*((1-t)*((1-t)*xOne + t*xTwo) + t*((1-t)*xTwo+t*xThree))
        y = (1-t)*((1-t)*((1-t)*yZero + t*yOne) + t*((1-t)*yOne+t*yTwo)) + t*((1-t)*((1-t)*yOne + t*yTwo) + t*((1-t)*yTwo+t*yThree))
        spline.append([x,y])
    return spline



def bezierIndex(positionVector,splineIndex,bezierSplinePoints):
    
    while positionVector[1] < bezierSplinePoints[splineIndex][1] and splineIndex < len(bezierSplinePoints)-1:
        
        splineIndex = splineIndex + 1

    return splineIndex




mainLoop()




fig = go.Figure()
fig.add_trace(go.Scatter(x=[point[0] for point in bezierSplinePoints],
                         y=[point[1] for point in bezierSplinePoints],
                         mode='markers',
                         marker=dict(color='blue', size=10),
                         ))

# Customize layout
fig.update_layout(title='Scatter Plot of Points',
                  xaxis=dict(title='X-axis' , range = [-500,500]),
                  yaxis=dict(title='Y-axis'),
                  )

# Show the plot
fig.show()









#with open(csvFileName, mode='w', newline='') as file:
#    writer = csv.writer(file)
#
#    writer.writerow(['Time Stamp', 'xPosition','yPosition','vehicleState'])
#
#    for i in range(len(logsTimeStamp)):
#
#        writer.writerow([logsTimeStamp[i], logsXPosition[i],logsYPosition[i],logsVehicleVector[i]])







timestamps = np.linspace(0, timeStep*1000, iterations)  # Example timestamps
x = logsXPosition
y = logsYPosition


magnifyFactor = 5
vehicleMagnifierFactor = 15
vehicleAirVeocityMagnifierFactor = .5

fig, ax = plt.subplots()
ax.set_aspect('equal')  # Set aspect ratio

text1 = ax.text(0.8, 0.9, '', transform=ax.transAxes, fontsize=12, ha='center')


def update(frame):

    normilizedVechicleVector = logsVehicleVector[frame]/np.linalg.norm(logsVehicleVector[frame])
    ax.set_title(f'Time: {timestamps[frame]:.2f}')  # Update title with timestamp
    vehicle = plt.arrow(x[frame], y[frame],vehicleMagnifierFactor*normilizedVechicleVector[0],-vehicleMagnifierFactor*normilizedVechicleVector[1], width = 2.5,edgecolor='none',facecolor ='green')

    lift = plt.arrow(x[frame], y[frame],logsLift[frame][0]*magnifyFactor,logsLift[frame][1]*magnifyFactor, width = 2,edgecolor='none',facecolor ='blue')

    drag = plt.arrow(x[frame], y[frame],logsDrag[frame][0]*magnifyFactor,logsDrag[frame][1]*magnifyFactor, width = 2,edgecolor='none', facecolor='red')

    finVector = plt.arrow(x[frame], y[frame],logsAppliedFinForce[frame][0]*4,logsAppliedFinForce[frame][1]*4, width = 1.5,edgecolor='none', facecolor='orange')

    #vehicleAirVelocity = plt.arrow(x[frame], y[frame],logsVehicleAirVelocity[frame][0]*vehicleAirVeocityMagnifierFactor,logsVehicleAirVelocity[frame][1]*vehicleAirVeocityMagnifierFactor, width = 2,edgecolor='none', facecolor='black')

    if frame %10 == 0:
        text1.set_text(f'Lift (blue) {logsLift[frame]}(N) \n Drag (red) {logsDrag[frame]}(N) \n Velocity {logsVelocity[frame]}(m/s)')
    return lift,drag,text1,vehicle,finVector,#vehicleAirVelocity,
# Create animation

interval = timeStep * 1000  # Interval between frames in milliseconds
animation = FuncAnimation(fig, update, frames=iterations, interval=interval, blit=True)

# Display the animation
plt.xlim(-300,300)  # Adjust limits as per your data range
plt.ylim(-10, 520)  # Adjust limits as per your data range
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Object Animation')

plt.show()
print(f"{logsRotation}\n")
#print(f"{logsTimeStamp}\n")
#print(f"{logsXPosition}\n")
#print(f"{logsYPosition}\n")


