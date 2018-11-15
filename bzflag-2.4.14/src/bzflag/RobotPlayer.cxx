/* bzflag
 * Copyright (c) 1993-2018 Tim Riker
 *
 * This package is free software;  you can redistribute it and/or
 * modify it under the terms of the license found in the file
 * named COPYING that should have accompanied this file.
 *
 * THIS PACKAGE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 *
 */

#include <array>
#include <map>
#include <vector>
#include "playing.h"

// interface header
#include "RobotPlayer.h"

// common implementation headers
#include "BZDBCache.h"

// local implementation headers
#include "World.h"
#include "Intersect.h"
#include "TargetingUtils.h"
#include "dectree.h"

std::vector<BzfRegion*>* RobotPlayer::obstacleList = NULL;

RobotPlayer::RobotPlayer(const PlayerId& _id, const char* _name,
			 ServerLink* _server,
			 const char* _motto = "") :
LocalPlayer(_id, _name, _motto),
target(NULL),
pathIndex(0),
timerForShot(0.0f),
drivingForward(true),
targetFlag(NULL),
flagsFound(0),
flocking(false),
seeker(false),
initialize_path(true),
pathfinder_mult(8),
need_path(true)

{
  gettingSound = false;
  server       = _server;
}

/* finds all colored flags in the world that are not the robots own color and stores pointer to them into enemyFlags vector.
 /*
 /* Added by Brendan McCabe
 */
void RobotPlayer::findEnemyFlags() {
  bool found = false;
  World *world = World::getWorld();
  for (int i = 0; i < world->getMaxFlags(); i++) {
    if (world->getFlag(i).type->flagTeam != getTeam() && world->getFlag(i).type->flagTeam != TeamColor::NoTeam) {
      enemyFlags.push_back(&(world->getFlag(i)));
      found = true;
    }
    if (found) {
      flagsFound = true;
    }
  }
}

/* description: finds all colored flags in the world that are not the robots own color and stores pointer to them into enemyFlags vector.
 /* params:      none
 /* return:      none
 /* Added by Brendan McCabe
 */
bool RobotPlayer::hasTargetFlag() {
  bool rValue = false;
  if (targetFlag != NULL) {
    if (targetFlag->type->getColor() == getFlag()->getColor()) {
      rValue = true;
    }
  }
  return rValue;
}

/* description: finds the distance between the robot tanks position and a flag
 /* params:      flag (Flag*): pointer to the flag who's distance the function finds
 /* return:      (float): the distance between the robot and flag
 /* Added by Brendan McCabe
 */
float RobotPlayer::flagDistance(const Flag *flag) {
  float xdist = getPosition()[0] - flag->position[0];
  float ydist = getPosition()[1] - flag->position[1];
  return (std::sqrt(std::pow(xdist, 2) + std::pow(ydist, 2)));
}

/* description: Sets the first element of enemyFlag as the robot's targeted flag
 /* params:      none
 /* return:      none
 /* Added by Brendan McCabe
 */
void RobotPlayer::seekFirstFlag() {
  targetFlag = enemyFlags.at(0);
}

/* description: Sets the enemy flag that's both the shortest distance away from
 /*                from the robot and not being carried by a tank as the robot's targeted flag.
 /* params:      none
 /* return:      none
 /* Added by Brendan McCabe
 */
void RobotPlayer::seekClosestFlag() {
  Flag* closestFlag = targetFlag;
  float closestDist = 9999999.9f;
  if (targetFlag == NULL) {
    targetFlag = enemyFlags.at(0);
  }
  if (getFlag() != targetFlag->type) {
    for (std::vector<Flag*>::iterator i = enemyFlags.begin(); i != enemyFlags.end(); i++) {
      float distance = flagDistance(*i);
      if (closestDist > distance && !((*i)->status == FlagStatus::FlagOnTank)) {
	closestDist = distance;
	closestFlag = *i;
      }
    }
  }
  targetFlag = closestFlag;
}

/* description: Sets the robot's targeted flag to be the flag of the highest scoring team whose f
 /*                flag is not being carried.
 /* params:      none
 /* return:      none
 /* Added by Brendan McCabe
 */
void RobotPlayer::seekLeaderFlag() {
  Flag* leaderFlag = targetFlag;
  World* world = World::getWorld();
  std::multimap<short, int> teamScores;
  for (int i = 0; i < NumTeams; i++) {
    Team* team = &(world->getTeam(i));
    if (team->size > 0) {
      teamScores.insert(std::make_pair(team->getWins() - team->getLosses(), i));
    }
  }
  bool searching = true;
  for (std::multimap<short, int>::reverse_iterator teamIter = teamScores.rbegin(); teamIter != teamScores.rend() && searching; teamIter++) {
    for (std::vector<Flag*>::reverse_iterator flagIter = enemyFlags.rbegin(); flagIter != enemyFlags.rend() && searching; flagIter++) {
      if ((*teamIter).second == (*flagIter)->type->flagTeam && (*flagIter)->status != FlagStatus::FlagOnTank) {
	leaderFlag = *flagIter;
	searching = false;
      }
    }
  }
  targetFlag = leaderFlag;
}

/* description: This function designates a robot as its team's seeker if there are no human players on
 /*              the team and it has the lowest player id on the team.  A seeker robot captures the flag rather
 /*              than flock
 /* params:      none
 /* return:      none
 /* Added by Brendan McCabe
 */
void RobotPlayer::findseeker() {
  bool seeker_found = false;
  std::vector<int> teamplayers;
  for (int t = 0; t <= World::getWorld()->getCurMaxPlayers(); t++)
  {
    Player *p = 0;
    if (t < World::getWorld()->getCurMaxPlayers())
      p = World::getWorld()->getPlayer(t);
    else
      p = LocalPlayer::getMyTank();
    if (p) {
      if (p->getTeam() == getTeam()) {
	teamplayers.push_back(t);
	if (!(p->getPlayerType() == PlayerType::ComputerPlayer)) {
	  seeker_found = true;
	}
      }
    }
  }
  if (!seeker_found && World::getWorld()->getPlayer(teamplayers.at(0))->getId() == getId()) {
    seeker = true;
  }
  else
    seeker = false;
}

/* description: computes the flocking vector for the flocking ai and then stores it in the object's flock direction
 /*              variable.  The flocking vector is calculated by adding the cohesion, seperation, and align vectors to
 /*              the current position.  This function contains constants which determine how significant each of the
 /*              the three vectors are and how close teammates must be for them to be applied
 /* params:      none
 /* return:      none
 /* Added by Brendan McCabe
 */
void RobotPlayer::flock() {
  const float COHESION_RANGE = BZDBCache::worldSize * 1.0f;
  const float SEPERATION_RANGE = BZDBCache::tankRadius * 13.0f;
  const float COHESION_MODIFIER = 12.0f;
  const float SEPERATION_MODIFIER = 9.0f;
  const float ALIGNMENT_MODIFIER = 17.0f;
  float cohesion[] = { 0.0, 0.0, 0.0 };
  float seperation[] = { 0.0, 0.0, 0.0 };
  float align[] = { 0.0f, 0.0f, 0.0f };
  flock_direction[0] = 0.0f;
  flock_direction[1] = 0.0f;
  calcCohesion(cohesion, COHESION_RANGE, COHESION_MODIFIER);
  calcSeperation(seperation, SEPERATION_RANGE, SEPERATION_MODIFIER);
  calcAlignment(align, COHESION_RANGE, ALIGNMENT_MODIFIER);
  flock_direction[0] = getPosition()[0] + cohesion[0] + seperation[0] + align[0];
  flock_direction[1] = getPosition()[1] + cohesion[1] + seperation[1] + align[0];
}

/* description: getter function for the targeted flag
 /* params:      none
 /* return:      (Flag*): pointer to the targeted flag
 /* Added by Brendan McCabe
 */
Flag* RobotPlayer::getTargetFlag() {
  return targetFlag;
}

/* description: Calculates the seperation vector of the flocking ai. Teammates who are closer
 /*              will have a greater affect on the vector.  When the modifier is 1.0, each teammate
 /*              will contribute a vector with a magnitude of the range minus the distance between the two tanks
 /* params:  seperation (float*): The array where the seperation vector will be stored
 /*          range (float): The maximum distance where the seperation vector will still be applied
 /*          modifier (float): multiplier for the strength of the vector
 /* return:      none
 /* Added by Brendan McCabe
 */
void RobotPlayer::calcSeperation(float* seperation, const float range, const float modifier) {
  TeamColor team = getTeam();
  int teammates = 0;
  seperation[0] = 0.0f;
  seperation[1] = 0.0f;
  for (int t = 0; t <= World::getWorld()->getCurMaxPlayers(); t++)
  {
    Player *p = 0;
    if (t == 0)
      p = LocalPlayer::getMyTank();
    else if (t < World::getWorld()->getCurMaxPlayers())
      p = World::getWorld()->getPlayer(t);
    else
      p = LocalPlayer::getMyTank();
    if (     p->getTeam() == team
	&& range > calcDistance(getPosition(), p->getPosition())
	&& getId() != p->getId()) {
      float angle = calcAngle(getPosition(), p->getPosition());
      teammates++;
      seperation[0] += (range * cos(angle)) - (getPosition()[0] - p->getPosition()[0]);
      seperation[1] += (range * sin(angle)) - (getPosition()[1] - p->getPosition()[1]);
    }
  }
  seperation[0] *= modifier;
  seperation[1] *= modifier;
}

/* description: Calculates the alignment vector of the flocking ai. The vector is the average velocity
 /*              of the team neighborhood
 /* params:  align (float*): The array where the alignment vector will be stored
 /*          range (float): The maximum distance where the alignment vector will still be applied
 /*          modifier (float): multiplier for the strength of the vector
 /* return:      none
 /* Added by Brendan McCabe
 */
void RobotPlayer::calcAlignment(float* align, const float range, const float modifier) {
  TeamColor team = getTeam();
  int teammates = 0;
  align[0] = 0.0f;
  align[1] = 0.0f;
  for (int t = 0; t <= World::getWorld()->getCurMaxPlayers(); t++)
  {
    Player *p = 0;
    if (t == 0)
      p = LocalPlayer::getMyTank();
    else if (t < World::getWorld()->getCurMaxPlayers())
      p = World::getWorld()->getPlayer(t);
    else
      p = LocalPlayer::getMyTank();
    if (p) {
      if (     p->getTeam() == team
	  && range > calcDistance(getPosition(), p->getPosition())) {
	teammates++;
	align[0] += p->getVelocity()[0];
	align[0] += p->getVelocity()[0];
      }
    }
  }
  align[0] /= (float)(teammates) * modifier;
  align[1] /= (float)(teammates) * modifier;
}

/* description: Calculates the cohesion vector of the flocking ai. Finds the team neiborhood's center of mass and then
 /*              subtracts the robots current position from it
 /* params:  align (float*): The array where the cohesion vector will be stored
 /*          range (float): The maximum distance where the cohesion vector will still be applied
 /*          modifier (float): multiplier for the strength of the vector
 /* return:      none
 /* Added by Brendan McCabe
 */
void RobotPlayer::calcCohesion(float* cohesion, const float range, const float modifier) {
  float com[] = { 0.0f, 0.0f, 0.0f };
  TeamColor team = getTeam();
  int teammates = 0;
  com[0] = 0.0f;
  com[1] = 0.0f;
  for (int t = 0; t <= World::getWorld()->getCurMaxPlayers(); t++)
  {
    Player *p = 0;
    if (t == 0)
      p = LocalPlayer::getMyTank();
    else if (t < World::getWorld()->getCurMaxPlayers())
      p = World::getWorld()->getPlayer(t);
    else
      p = LocalPlayer::getMyTank();
    if (p) {
      if (     p->getTeam() == team
	  && range > calcDistance(getPosition(), p->getPosition())) {
	teammates++;
	com[0] += p->getPosition()[0];
	com[1] += p->getPosition()[1];
      }
    }
  }
  com[0] /= (float)(teammates);
  com[1] /= (float)(teammates);
  cohesion[0] = modifier * (com[0] - getPosition()[0]);
  cohesion[1] = modifier * (com[1] - getPosition()[1]);
}

/* description: Calculates the distance between two two component vectors
 /* params:      a (float*) pointer to the first vector
 b (float*) pointer to the second vector
 /* return:      (float) the distance
 /* Added by Brendan McCabe
 */
float RobotPlayer::calcDistance(const float* a, const float* b) {
  return (sqrt(pow((a[0] - b[0]), 2) + pow((a[1] - b[1]), 2)));
}
float RobotPlayer::calcDistance(float* a, const float* b) {
  return (sqrt(pow((a[0] - b[0]), 2) + pow((a[1] - b[1]), 2)));
}
float RobotPlayer::calcDistance(const int* a, const float* b) {
  return (sqrt(pow(((float)a[0] - b[0]), 2) + pow(((float)a[1] - b[1]), 2)));
}

/* description: Calculates the angle between two two component vectors
 /* params:      a (float*) pointer to the first vector
 b (float*) pointer to the second vector
 /* return:      (float) the angle
 /* Added by Brendan McCabe
 */
float RobotPlayer::calcAngle(const float* a, const float* b) {
  return atan2((a[1] - b[1]), (a[0] - b[0]));
}


/*
 Function name:	initialize_pathfinder
 Parameters:     none
 Return:		none
 Description:	allocates memory for the pathfinder object
 /* Added by Brendan McCabe
 */
void RobotPlayer::initialize_pathfinder() {
  pathfinder.set_dimensions((int)(BZDBCache::worldSize / 2), (int)(BZDBCache::worldSize / 2));
  pathfinder.build_graph();
  for (int i = 0; i < (BZDBCache::worldSize * 2 / pathfinder_mult); i++) {
    for (int j = 0; j < (BZDBCache::worldSize * 2 / pathfinder_mult); j++) {
      float pos[3];
      pos[0] = pathfinder_mult*i - BZDBCache::worldSize;
      pos[1] = pathfinder_mult*j - BZDBCache::worldSize;
      pos[2] = 0;
      if (World::getWorld()->inBuilding(pos, BZDBCache::tankRadius / 2, BZDBCache::tankHeight) != NULL) {
	pathfinder.set_obstacle(i, j);
      }
    }
  }
  aicore::DecisionTrees::init();
}

/*
 Function name:	find_path
 Parameters:     none
 Return:		none
 Description:	finds the shortest path to the targeted flag or base through the pathfinder object
 path is stored in the my_path variable
 /* Added by Brendan McCabe
 */
void RobotPlayer::find_path() {
  need_path = false;
  pathfinder.clear_graph();
  pathfinder.set_location((int)((getPosition()[0] + BZDBCache::worldSize) / pathfinder_mult), (int)((getPosition()[1] + BZDBCache::worldSize) / pathfinder_mult));
  if (hasTargetFlag()) {
    const float *target_pos = World::getWorld()->getBase(getTeam(), 0);
    pathfinder.set_target((int)((target_pos[0] + BZDBCache::worldSize) / (int)(pathfinder_mult)), (int)((target_pos[1] + BZDBCache::worldSize))/ (int)(pathfinder_mult));
  }
  else {
    const float *target_pos = targetFlag->position;
    pathfinder.set_target((int)((target_pos[0] + BZDBCache::worldSize) / (int)(pathfinder_mult)), (int)((target_pos[1] + BZDBCache::worldSize)) / (int)(pathfinder_mult));
  }
  //controlPanel->addMessage(me);
  pathfinder.set_distances(10, 14);
  pathfinder.find_path();
  my_path.clear();
  MyPathfinder::Cell* path_walker = pathfinder.target;
  while (path_walker != NULL){
    std::array<int, 3> next_cell;
    next_cell[0] = path_walker->position[0];
    next_cell[1] = path_walker->position[1];
    next_cell[2] = 0;
    my_path.insert(my_path.begin(), next_cell);
    path_walker = path_walker->parent;
  }
}

/*
 Function name:	follow_path
 Parameters:     none
 Return:		none
 Description:	stores location of the next cell in my_path in the destination variable which
 will be used to tell the robot where to go
 /* Added by Brendan McCabe
 */
void RobotPlayer::follow_path() {
  float next_coor[3];
  if (!my_path.empty()) {
    next_coor[0] = (my_path.at(0).data()[0] * pathfinder_mult) - BZDBCache::worldSize;
    next_coor[1] = (my_path.at(0).data()[1] * pathfinder_mult) - BZDBCache::worldSize;
    next_coor[2] = (my_path.at(0).data()[2] * pathfinder_mult) - BZDBCache::worldSize;
  }
  if (calcDistance(next_coor, getPosition()) < BZDBCache::tankRadius * 3 && my_path.size() > 1) {
    my_path.erase(my_path.begin());
  }
  if (my_path.size() == 1) {
    if (has_reached()) {
      my_path.erase(my_path.begin());
    }
    if (hasTargetFlag()) {
      destination[0] = World::getWorld()->getBase(getTeam(), 0)[0];
      destination[1] = World::getWorld()->getBase(getTeam(), 0)[1];
    }
    else {
      destination[0] = targetFlag->position[0];
      destination[1] = targetFlag->position[1];
    }
  }
  else if (!my_path.empty()) {
    destination[0] = (float)(my_path.at(0).data()[0] * pathfinder_mult) - BZDBCache::worldSize;
    destination[1] = (float)(my_path.at(0).data()[1] * pathfinder_mult) - BZDBCache::worldSize;
  }
  else {
    need_path = true;
    destination[0] = getPosition()[0];
    destination[1] = getPosition()[1];
  }
  destination[2] = 0;
}

/*
 Function name:	has_reached
 Parameters:
 Return:		true: tank has reached its objective or is far enough away that it needs to find a new path
 false: tank is close to objective but has not reached it
 Description:	used to determine if the tank has reached the objective with one cell remaining on the path
 /* Added by Brendan McCabe
 */
bool RobotPlayer::has_reached() {
  if (calcDistance(targetFlag->position, getPosition()) < BZDBCache::tankRadius * pathfinder_mult && !hasTargetFlag()) {
    return false;
  }
  else if (calcDistance(World::getWorld()->getBase(getTeam(), 0), getPosition()) < BZDBCache::tankRadius * pathfinder_mult && hasTargetFlag()) {
    return false;
  }
  return true;
}

// estimate a player's position at now+t, similar to dead reckoning
void RobotPlayer::projectPosition(const Player *targ,const float t,float &x,float &y,float &z) const
{
  double hisx=targ->getPosition()[0];
  double hisy=targ->getPosition()[1];
  double hisz=targ->getPosition()[2];
  double hisvx=targ->getVelocity()[0];
  double hisvy=targ->getVelocity()[1];
  double hisvz=targ->getVelocity()[2];
  double omega=fabs(targ->getAngularVelocity());
  double sx,sy;
  
  if ((targ->getStatus() & PlayerState::Falling) || fabs(omega) < 2*M_PI / 360 * 0.5)
  {
    sx=t*hisvx;
    sy=t*hisvy;
  }
  else
  {
    double hisspeed = hypotf(hisvx, hisvy);
    double alfa = omega * t;
    double r = hisspeed / fabs(omega);
    double dx = r * sin(alfa);
    double dy2 = r * (1 - cos(alfa));
    double beta = atan2(dy2, dx) * (targ->getAngularVelocity() > 0 ? 1 : -1);
    double gamma = atan2(hisvy, hisvx);
    double rho = gamma+beta;
    sx = hisspeed * t * cos(rho);
    sy = hisspeed * t * sin(rho);
  }
  x=(float)hisx+(float)sx;
  y=(float)hisy+(float)sy;
  z=(float)hisz+(float)hisvz*t;
  if (targ->getStatus() & PlayerState::Falling)
    z += 0.5f * BZDBCache::gravity * t * t;
  if (z < 0) z = 0;
}

/*****start Decision tree action and decision functions added by aidan akamine *******/
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*
 * same as Player::isAlive(), needed to match type of decision tree
 */
bool        RobotPlayer::amAlive(float dt)
{
  return isAlive();
}

void        RobotPlayer::doNothing(float dt)
{
  
}


/******************************* is robot about to get hit *******************************/
/******************************* is robot about to get hit ***********************************/
bool            RobotPlayer::isInDanger(float dt)
{
  
  for (int t=0; t <= World::getWorld()->getCurMaxPlayers(); t++)
  {
    Player *p = 0;
    if (t < World::getWorld()->getCurMaxPlayers())
      p = World::getWorld()->getPlayer(t);
    else
      p = LocalPlayer::getMyTank();
    if (!p || p->getId() == getId())
      continue;
    const int maxShots = p->getMaxShots();
    for (int s = 0; s < maxShots; s++) {
      ShotPath* shot = p->getShot(s);
      if (!shot || shot->isExpired())
	continue;
      // ignore invisible bullets completely for now (even when visible)
      if (shot->getFlag() == Flags::InvisibleBullet)
	return false;
      const float* position = getPosition();
      shotPos = shot->getPosition();
      if ((fabs(shotPos[2] - position[2]) > BZDBCache::tankHeight) && (shot->getFlag() != Flags::GuidedMissile))
	continue;
      const float dist = TargetingUtils::getTargetDistance(position, shotPos);
      if (dist < 150.0f) {
	shotVel = shot->getVelocity();
	float shotAngle = atan2f(shotVel[1], shotVel[0]);
	float shotUnitVec[2] = {cosf(shotAngle), sinf(shotAngle)};
	
	float trueVec[2] = {(position[0]-shotPos[0])/dist,(position[1]-shotPos[1])/dist};
	float dotProd = trueVec[0]*shotUnitVec[0]+trueVec[1]*shotUnitVec[1];
	if (dotProd > 0.97f) {
	  evading = true;
	  return true;
	}
      }
    }
  }
  evading = false;
  return false;
}

/******************************* check firing status**********************************************/
/******************************* check firing status**********************************************/
bool            RobotPlayer::readyToFire(float dt)
{
  timerForShot  -= dt;
  if (timerForShot < 0.0f){
    timerForShot = 0.0f;
  }
  if (getFiringStatus() == Ready) {
    return true;
  }
  return false;
}


/******************************* check shot timer**********************************************/
/******************************* check shot timer**********************************************/
bool            RobotPlayer::shotTimerReady(float dt)
{
  if (timerForShot <= 0.0f){
    timerForShot = 0.0f;
    if(!path.empty()){
      return true;
    }
  }
  return false;
}

/******************************* check is shot is within half a tank length*******************************************/
/******************************* check is shot is within half a tank length*******************************************/
bool            RobotPlayer::shotWithinHalfTankLenght(float dt)
{
  
  float tankRadius = BZDBCache::tankRadius;
  const float shotRadius = BZDB.eval(StateDatabase::BZDB_SHOTRADIUS);
  float p1[3];
  const float azimuth = getAngle();
  getProjectedPosition(target, p1);
  const float* p2 = getPosition();
  shootingAngle = atan2f(p1[1] - p2[1], p1[0] - p2[0]);
  if (shootingAngle < 0.0f)
    shootingAngle += (float)(2.0 * M_PI);
  float azimuthDiff   = shootingAngle - azimuth;
  if (azimuthDiff > M_PI)
    azimuthDiff -= (float)(2.0 * M_PI);
  else
    if (azimuthDiff < -M_PI)
      azimuthDiff += (float)(2.0 * M_PI);
  
  const float targetdistance = hypotf(p1[0] - p2[0], p1[1] - p2[1]) -
  BZDB.eval(StateDatabase::BZDB_MUZZLEFRONT) - tankRadius;
  
  const float missby = fabs(azimuthDiff) * (targetdistance - BZDBCache::tankLength);
  // only shoot if we miss by less than half a tanklength and no building inbetween
  if (missby < 0.5f * BZDBCache::tankLength && p1[2] < shotRadius) {
    return true;
  }
  return false;
}



/******************************* check if building is between shot *******************************************/
/******************************* check if building is between shot *******************************************/
bool            RobotPlayer::isBuildingInBetween(float dt)
{
  
  float tankRadius = BZDBCache::tankRadius;
  float p1[3];
  const float* p2 = getPosition();
  const float azimuth = getAngle();
  getProjectedPosition(target, p1);
  float pos[3] = {getPosition()[0], getPosition()[1], getPosition()[2] +  BZDB.eval(StateDatabase::BZDB_MUZZLEHEIGHT)};
  float dir[3] = {cosf(azimuth), sinf(azimuth), 0.0f};
  Ray tankRay(pos, dir);
  const float targetdistance = hypotf(p1[0] - p2[0], p1[1] - p2[1]) - BZDB.eval(StateDatabase::BZDB_MUZZLEFRONT) - tankRadius;
  float maxdistance = targetdistance;
  if (!ShotStrategy::getFirstBuilding(tankRay, -0.5f, maxdistance))
  {
    return true;
  }
  return false;
}

/******************************* check is shot is will hit a teammate *******************************************/
/******************************* check is shot is will hit a teammate *******************************************/
bool            RobotPlayer::WillShotHitTeammate(float dt)
{
  
  const float shotRange  = BZDB.eval(StateDatabase::BZDB_SHOTRANGE);
  const float azimuth = getAngle();
  float pos[3] = {getPosition()[0], getPosition()[1], getPosition()[2] +  BZDB.eval(StateDatabase::BZDB_MUZZLEHEIGHT)};
  float dir[3] = {cosf(azimuth), sinf(azimuth), 0.0f};
  
  Ray tankRay(pos, dir);
  // try to not aim at teammates
  for (int i=0; i <= World::getWorld()->getCurMaxPlayers(); i++)
  {
    Player *p = 0;
    if (i < World::getWorld()->getCurMaxPlayers())
      p = World::getWorld()->getPlayer(i);
    else
      p = LocalPlayer::getMyTank();
    
    if (!p || p->getId() == getId() || validTeamTarget(p) ||
	!p->isAlive()) continue;
    float relpos[3] = {getPosition()[0] - p->getPosition()[0],
      getPosition()[1] - p->getPosition()[1],
      getPosition()[2] - p->getPosition()[2]};
    Ray ray(relpos, dir);
    float impact = rayAtDistanceFromOrigin(ray, 5 * BZDBCache::tankRadius);
    if (impact > 0 && impact < shotRange)
    {
      return false;
    }
  }
  
  return true;
}



/******************************* check if tank is holding a flag *******************************************/
/******************************* check if tank is holding a flag *******************************************/
bool             RobotPlayer::holdingFlag(float dt)
{
  if(this->getFlag() != Flags::Null){
    return true;
  }
  return false;
}



/******************************* check if tank is holding a sticky flag *******************************************/
/******************************* check if tank is holding a sticky flag *******************************************/
bool             RobotPlayer::isFlagSticky(float dt)
{
  if(this->getFlag()->endurance == FlagSticky){
    return true;
  }
  return false;
}


/******************************* check if tank is holding a team flag *******************************************/
/******************************* check if tank is holding a team flag *******************************************/
bool             RobotPlayer::isFlagTeamFlag(float dt)
{
  if(this->getFlag()->flagTeam != NoTeam){
    return true;
  }
  return false;
}


/******************************* check if tank is holding their own flag *******************************************/
/******************************* check if tank is holding their own flag *******************************************/
bool             RobotPlayer::holdingMyTeamFlag(float dt)
{
  if(this->getFlag()->flagTeam == getTeam()){
    return true;
  }
  return false;
}



/******************************* drop flag *******************************************/
/******************************* drop flag *******************************************/
void             RobotPlayer::dropFlag(float dt)
{
  serverLink->sendDropFlag(getId(), getPosition());
}

/*
 if(this->getFlag()->flagTeam == getTeam()){
 serverLink->sendDropFlag(getId(), position);
 } else if (this->getFlag()->flagTeam == TeamColor::NoTeam){
 if(this->getFlag()->endurance != FlagSticky){
 serverLink->sendDropFlag(getId(), position);
 }
 }*/

/******************************* tank is firing shot *******************************************/
/******************************* tank is firing shot *******************************************/
void             RobotPlayer::Shooting(float dt)
{
  if (fireShot()) {
    timerForShot = float(bzfrand()) * 0.6f + 0.2f;
  }
}


/******************************* teammate is blocking shot *******************************************/
/******************************* teammate is blocking shot *******************************************/
void             RobotPlayer::teammateBlockingShot(float dt)
{
  timerForShot = 0.1f;
}

/******************************* have robot evade    **********************************************/
/******************************* have robot evade    **********************************************/
void            RobotPlayer::evade(float dt)
{
  /*
   char buffer[280];
   sprintf (buffer, "evade" );
   controlPanel->addMessage("evading");
   */
  const float azimuth = getAngle();
  const float* position = getPosition();
  const float dist = TargetingUtils::getTargetDistance(position, shotPos);
  float shotAngle = atan2f(shotVel[1], shotVel[0]);
  float shotUnitVec[2] = {cosf(shotAngle), sinf(shotAngle)};
  
  float trueVec[2] = {(position[0]-shotPos[0])/dist,(position[1]-shotPos[1])/dist};
  float dotProd = trueVec[0]*shotUnitVec[0]+trueVec[1]*shotUnitVec[1];
  
  if (dotProd > 0.97f) {
    float rotation;
    float rotation1 = (float)((shotAngle + M_PI/2.0) - azimuth);
    if (rotation1 < -1.0f * M_PI) rotation1 += (float)(2.0 * M_PI);
    if (rotation1 > 1.0f * M_PI) rotation1 -= (float)(2.0 * M_PI);
    
    float rotation2 = (float)((shotAngle - M_PI/2.0) - azimuth);
    if (rotation2 < -1.0f * M_PI) rotation2 += (float)(2.0 * M_PI);
    if (rotation2 > 1.0f * M_PI) rotation2 -= (float)(2.0 * M_PI);
    
    if (fabs(rotation1) < fabs(rotation2))
      rotation = rotation1;
    else
      rotation = rotation2;
    setDesiredSpeed(1.0f);
    setDesiredAngVel(rotation);
    evading = true;
  }
  
}


void            RobotPlayer::followPath(float dt)
{
  
  float epInit[] = { 0.0, 0.0, 0.0 };
  endPoint = epInit;
  //findseeker();
  if (!flagsFound) {
    findEnemyFlags();
  }
  seekClosestFlag();
  
  if (need_path) {
    find_path();
  }
  follow_path();
  endPoint = destination;
  /*End of code modified by Brendan McCabe*/
  
}
/*****end Decision tree action and decision functions added by aidan akamine *******/
/*******************************************************************************/


// get coordinates to aim at when shooting a player; steps:
// 1. estimate how long it will take shot to hit target
// 2. calc position of target at that point of time
// 3. jump to 1., using projected position, loop until result is stable
void RobotPlayer::getProjectedPosition(const Player *targ, float *projpos) const
{
  double myx = getPosition()[0];
  double myy = getPosition()[1];
  double hisx = targ->getPosition()[0];
  double hisy = targ->getPosition()[1];
  double deltax = hisx - myx;
  double deltay = hisy - myy;
  double distance = hypotf(deltax,deltay) - BZDB.eval(StateDatabase::BZDB_MUZZLEFRONT) - BZDBCache::tankRadius;
  if (distance <= 0) distance = 0;
  double shotspeed = BZDB.eval(StateDatabase::BZDB_SHOTSPEED)*
  (getFlag() == Flags::Laser ? BZDB.eval(StateDatabase::BZDB_LASERADVEL) :
   getFlag() == Flags::RapidFire ? BZDB.eval(StateDatabase::BZDB_RFIREADVEL) :
   getFlag() == Flags::MachineGun ? BZDB.eval(StateDatabase::BZDB_MGUNADVEL) : 1) +
  hypotf(getVelocity()[0], getVelocity()[1]);
  
  double errdistance = 1.0;
  float tx, ty, tz;
  for (int tries=0 ; errdistance > 0.05 && tries < 4 ; tries++)
  {
    float t = (float)distance / (float)shotspeed;
    projectPosition(targ, t + 0.05f, tx, ty, tz); // add 50ms for lag
    double distance2 = hypotf(tx - myx, ty - myy);
    errdistance = fabs(distance2-distance) / (distance + ZERO_TOLERANCE);
    distance = distance2;
  }
  projpos[0] = tx; projpos[1] = ty; projpos[2] = tz;
  
  // projected pos in building -> use current pos
  if (World::getWorld()->inBuilding(projpos, 0.0f, BZDBCache::tankHeight)) {
    projpos[0] = targ->getPosition()[0];
    projpos[1] = targ->getPosition()[1];
    projpos[2] = targ->getPosition()[2];
  }
}


void			RobotPlayer::doUpdate(float dt)
{
  LocalPlayer::doUpdate(dt);
  if (initialize_path) {
    initialize_path = false;
    initialize_pathfinder();
  }
  /*lines modified by aidan akamine*/
  
  aicore::DecisionTreeNode *node1 = aicore::DecisionTrees::shootingDecisions[0].makeDecision(this, dt);
  (this->*(((aicore::ActionPtr*)node1)->actFuncPtr))(dt);
  aicore::DecisionTreeNode *node2 = aicore::DecisionTrees::holdingFlagDecisions[0].makeDecision(this, dt);
  (this->*(((aicore::ActionPtr*)node2)->actFuncPtr))(dt);
  
  /*end of lines modified by aidan akamine*/
  
  float tankRadius = BZDBCache::tankRadius;
  const float shotRange  = BZDB.eval(StateDatabase::BZDB_SHOTRANGE);
  const float shotRadius = BZDB.eval(StateDatabase::BZDB_SHOTRADIUS);
}

void			RobotPlayer::doUpdateMotion(float dt)
{
  if (initialize_path) {
    initialize_path = false;
    initialize_pathfinder();
  }
  aicore::DecisionTreeNode *node = aicore::DecisionTrees::doUpdateMotionDecisions[0].makeDecision(this, dt);
  (this->*(((aicore::ActionPtr*)node)->actFuncPtr))(dt);
  
  /*The following code was modified by Brendan McCabe*/
  if (isAlive()) {
    
    // record previous position
    const float oldAzimuth = getAngle();
    const float* oldPosition = getPosition();
    float position[3];
    position[0] = oldPosition[0];
    position[1] = oldPosition[1];
    position[2] = oldPosition[2];
    float azimuth = oldAzimuth;
    float tankAngVel = BZDB.eval(StateDatabase::BZDB_TANKANGVEL);
    float tankSpeed = BZDBCache::tankSpeed;
    
    // when we are not evading, follow the path
    if (!evading && dt > 0.0 && pathIndex < (int)path.size()) {
      float distance;
      float v[2];
      //endPoint = path[pathIndex].get();
      
      // find how long it will take to get to next path segment
      v[0] = endPoint[0] - position[0];
      v[1] = endPoint[1] - position[1];
      distance = hypotf(v[0], v[1]);
      float tankRadius = BZDBCache::tankRadius;
      // smooth path a little by turning early at corners, might get us stuck, though
      if (distance <= 2.5f * tankRadius)
	pathIndex++;
      
      float segmentAzimuth = atan2f(v[1], v[0]);
      float azimuthDiff = segmentAzimuth - azimuth;
      if (azimuthDiff > M_PI) azimuthDiff -= (float)(2.0 * M_PI);
      else if (azimuthDiff < -M_PI) azimuthDiff += (float)(2.0 * M_PI);
      if (fabs(azimuthDiff) > 0.01f) {
	// drive backward when target is behind, try to stick to last direction
	if (drivingForward)
	  drivingForward = fabs(azimuthDiff) < M_PI/2*0.9 ? true : false;
	else
	  drivingForward = fabs(azimuthDiff) < M_PI/2*0.3 ? true : false;
	setDesiredSpeed(drivingForward ? 1.0f : -1.0f);
	// set desired turn speed
	if (azimuthDiff >= dt * tankAngVel) {
	  setDesiredAngVel(1.0f);
	} else if (azimuthDiff <= -dt * tankAngVel) {
	  setDesiredAngVel(-1.0f);
	} else {
	  setDesiredAngVel(azimuthDiff / dt / tankAngVel);
	}
      } else {
	drivingForward = true;
	// tank doesn't turn while moving forward
	setDesiredAngVel(0.0f);
	// find how long it will take to get to next path segment
	if (distance <= dt * tankSpeed) {
	  pathIndex++;
	  // set desired speed
	  setDesiredSpeed(distance / dt / tankSpeed);
	}
	else if (seeker) {          // modified by
	  setDesiredSpeed(0.85f);   //Brendan McCabe
	}
	else  {
	  setDesiredSpeed(1.0f);
	}
      }
    }
  }
  LocalPlayer::doUpdateMotion(dt);
}

void			RobotPlayer::explodeTank()
{
  LocalPlayer::explodeTank();
  target = NULL;
  need_path = true;
  path.clear();
}

void			RobotPlayer::restart(const float* pos, float _azimuth)
{
  LocalPlayer::restart(pos, _azimuth);
  // no target
  path.clear();
  target = NULL;
  need_path = true;
  pathIndex = 0;
  
}

float			RobotPlayer::getTargetPriority(const
						       Player* _target) const
{
  // don't target teammates or myself
  if (!this->validTeamTarget(_target))
    return 0.0f;
  
  // go after closest player
  // FIXME -- this is a pretty stupid heuristic
  const float worldSize = BZDBCache::worldSize;
  const float* p1 = getPosition();
  const float* p2 = _target->getPosition();
  
  float basePriority = 1.0f;
  // give bonus to non-paused player
  if (!_target->isPaused())
    basePriority += 2.0f;
  // give bonus to non-deadzone targets
  if (obstacleList) {
    float nearest[2];
    const BzfRegion* targetRegion = findRegion (p2, nearest);
    if (targetRegion && targetRegion->isInside(p2))
      basePriority += 1.0f;
  }
  return basePriority
  - 0.5f * hypotf(p2[0] - p1[0], p2[1] - p1[1]) / worldSize;
}

void		    RobotPlayer::setObstacleList(std::vector<BzfRegion*>*
						 _obstacleList)
{
  obstacleList = _obstacleList;
}

const Player*		RobotPlayer::getTarget() const
{
  return target;
}

void			RobotPlayer::setTarget(const Player* _target)
{
  static int mailbox = 0;
  
  path.clear();
  target = _target;
  if (!target) return;
  
  // work backwards (from target to me)
  float proj[3];
  getProjectedPosition(target, proj);
  const float *p1 = proj;
  const float* p2 = getPosition();
  float q1[2], q2[2];
  BzfRegion* headRegion = findRegion(p1, q1);
  BzfRegion* tailRegion = findRegion(p2, q2);
  if (!headRegion || !tailRegion) {
    // if can't reach target then forget it
    return;
  }
  
  mailbox++;
  headRegion->setPathStuff(0.0f, NULL, q1, mailbox);
  RegionPriorityQueue queue;
  queue.insert(headRegion, 0.0f);
  BzfRegion* next;
  while (!queue.isEmpty() && (next = queue.remove()) != tailRegion)
    findPath(queue, next, tailRegion, q2, mailbox);
  
  // get list of points to go through to reach the target
  next = tailRegion;
  do {
    p1 = next->getA();
    path.push_back(p1);
    next = next->getTarget();
  } while (next && next != headRegion);
  if (next || tailRegion == headRegion)
    path.push_back(q1);
  else
    path.clear();
  pathIndex = 0;
}

BzfRegion*		RobotPlayer::findRegion(const float p[2],
						float nearest[2]) const
{
  nearest[0] = p[0];
  nearest[1] = p[1];
  const int count = obstacleList->size();
  for (int o = 0; o < count; o++)
    if ((*obstacleList)[o]->isInside(p))
      return (*obstacleList)[o];
  
  // point is outside: find nearest region
  float      distance      = maxDistance;
  BzfRegion* nearestRegion = NULL;
  for (int i = 0; i < count; i++) {
    float currNearest[2];
    float currDistance = (*obstacleList)[i]->getDistance(p, currNearest);
    if (currDistance < distance) {
      nearestRegion = (*obstacleList)[i];
      distance = currDistance;
      nearest[0] = currNearest[0];
      nearest[1] = currNearest[1];
    }
  }
  return nearestRegion;
}

float			RobotPlayer::getRegionExitPoint(
							const float p1[2], const float p2[2],
							const float a[2], const float targetPoint[2],
							float mid[2], float& priority)
{
  float b[2];
  b[0] = targetPoint[0] - a[0];
  b[1] = targetPoint[1] - a[1];
  float d[2];
  d[0] = p2[0] - p1[0];
  d[1] = p2[1] - p1[1];
  
  float vect = d[0] * b[1] - d[1] * b[0];
  float t    = 0.0f;  // safe value
  if (fabs(vect) > ZERO_TOLERANCE) {
    // compute intersection along (p1,d) with (a,b)
    t = (a[0] * b[1] - a[1] * b[0] - p1[0] * b[1] + p1[1] * b[0]) / vect;
    if (t > 1.0f)
      t = 1.0f;
    else if (t < 0.0f)
      t = 0.0f;
  }
  
  mid[0] = p1[0] + t * d[0];
  mid[1] = p1[1] + t * d[1];
  
  const float distance = hypotf(a[0] - mid[0], a[1] - mid[1]);
  // priority is to minimize distance traveled and distance left to go
  priority = distance + hypotf(targetPoint[0] - mid[0], targetPoint[1] - mid[1]);
  // return distance traveled
  return distance;
}

void			RobotPlayer::findPath(RegionPriorityQueue& queue,
					      BzfRegion* region,
					      BzfRegion* targetRegion,
					      const float targetPoint[2],
					      int mailbox)
{
  const int numEdges = region->getNumSides();
  for (int i = 0; i < numEdges; i++) {
    BzfRegion* neighbor = region->getNeighbor(i);
    if (!neighbor) continue;
    
    const float* p1 = region->getCorner(i).get();
    const float* p2 = region->getCorner((i+1)%numEdges).get();
    float mid[2], priority;
    float total = getRegionExitPoint(p1, p2, region->getA(),
				     targetPoint, mid, priority);
    priority += region->getDistance();
    if (neighbor == targetRegion)
      total += hypotf(targetPoint[0] - mid[0], targetPoint[1] - mid[1]);
    total += region->getDistance();
    if (neighbor->test(mailbox) || total < neighbor->getDistance()) {
      neighbor->setPathStuff(total, region, mid, mailbox);
      queue.insert(neighbor, priority);
    }
  }
}


// Local Variables: ***
// mode: C++ ***
// tab-width: 8 ***
// c-basic-offset: 2 ***
// indent-tabs-mode: t ***
// End: ***
// ex: shiftwidth=2 tabstop=8
