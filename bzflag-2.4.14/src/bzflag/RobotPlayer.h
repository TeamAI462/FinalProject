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

#ifndef	BZF_ROBOT_PLAYER_H
#define	BZF_ROBOT_PLAYER_H

#include "common.h"

/* system interface headers */
#include <vector>
#include <array>

/* interface header */
#include "LocalPlayer.h"

/* local interface headers */
#include "Region.h"
#include "RegionPriorityQueue.h"
#include "ServerLink.h"

#include "MyPathfinder.h"
#include <Stack>


class RobotPlayer : public LocalPlayer {
  public:
  RobotPlayer(const PlayerId&,
	      const char* name, ServerLink*,
	      const char* _motto);
  
  /*start of functions writen by Brendan McCabe*/
  
  //assignment 1 functions
  void                findEnemyFlags();
  bool                hasTargetFlag();
  void                seekFirstFlag();
  void                seekClosestFlag();
  void                seekLeaderFlag();
  void 		      seekWeakestFlag();
  Flag*               getTargetFlag();
  void                flock();
  std::vector<Flag*>  enemyFlags;
  static const int    defendingtanks = 2;
  static const int    Attackingtanks = 3;
  static int 	      numAttacking;
  static int	      numDefending;
  static const int    samePositionTimerMax = 400;
  
  //assignment 2 functions
  void                calcCohesion(float* cohesion, const float range, const float modifier);
  float               calcDistance(const float* a, const float* b);
  float               calcDistance(float* a, const float* b);
  float               calcDistance(const int* a, const float* b);
  float               calcAngle(const float* a, const float* b);
  void                calcSeperation(float* seperation, const float range, const float modifier);
  void                calcAlignment(float* align, const float range, const float modifier);
  void                findseeker();
  
  //assignment 3 functions
  void                initialize_pathfinder();
  void                find_path();
  void                follow_path();
/*added by aidan akamine*/
  void		      follow_path_to_position();
  void		      follow_path_to_unstuck();
  void		      find_Unstuck_path();
  void		      find_position_path();
  void		      setUnstuckDestination();
  void 		      setPositionToBase();
  void		      followPathToUnstuck();
  void		      becomeUnstuck(float dt);
/*end of lines added by aidan akamine*/
  bool                has_reached();
  
/*start of decision trees written by aidan akamine
  ///************************decision trees************************/
  void		doNothing(float dt);
  bool		amAlive(float dt);
  void 		assignPosition();
  /* movement decision tree***********************/
  /*choices*/
  bool		isInDanger(float dt);
  bool		isAttacking(float dt);
  void		evade(float dt);
  const float *	shotVel;
  const float *	shotPos;
  bool		evading;
  bool		stuck;
  int 		unstuck_steps;
  bool		goToBase;
  bool		goToMiddle;
  bool 		isStuck(float dt);
  bool		isNearBase(float dt);
  void		attackFlags(float dt);
  void		defendFlag(float dt);
  
  /*actions*/
  void 		followPath(float dt);
  void		goToPosition(float dt);
  
  /* shooting decision tree***********************/
  /*choices*/
  bool		readyToFire(float dt);
  bool		shotTimerReady(float dt);
  bool		shotWithinHalfTankLenght(float dt);
  float		shootingAngle;
  bool		isBuildingInBetween(float dt);
  bool		WillShotHitTeammate(float dt);
  /*actions*/
  void		Shooting(float dt);
  void 		teammateBlockingShot(float dt);
  
  /*dropping flag decision tree ******************/
  /*choices*/
  bool		holdingFlag(float dt);
  bool		isFlagSticky(float dt);
  bool		isFlagTeamFlag(float dt);
  bool		holdingMyTeamFlag(float dt);
  bool		holdingGoodFlag(float dt);
  bool 		targetFlagWithinRange(float dt);
  /*actions*/
  void		dropFlag(float dt);
/*end of decision trees written by aidan akamine
  
/*start of decision tree written by cassandra largosa
  /* seeking flag decision tree ******************/
  /**** RobotPlayer.h flag targeting decision tree ****/
  // choices
  bool          winning(float dt);
  bool          equalLoses(float dt);
  // actions
  void          targetClosest(float dt);
  void          targetLeader(float dt);
  void          targetWeakest(float dt);
/* end of decision tree written by cassandra largosa
  
  /****************end of functions written by Brendan McCabe**********************/
  
  float		getTargetPriority(const Player*) const;
  const Player*	getTarget() const;
  void		setTarget(const Player*);
  static void		setObstacleList(std::vector<BzfRegion*>*);
  
  void		restart(const float* pos, float azimuth);
  void		explodeTank();
  
  private:
  float               flagDistance(const Flag *flag); //added by Brendan McCabe
  void		doUpdate(float dt);
  void		doUpdateMotion(float dt);
  BzfRegion*		findRegion(const float p[2], float nearest[2]) const;
  float		getRegionExitPoint(
				   const float p1[2], const float p2[2],
				   const float a[2], const float targetPoint[2],
				   float mid[2], float& priority);
  void		findPath(RegionPriorityQueue& queue,
			 BzfRegion* region, BzfRegion* targetRegion,
			 const float targetPoint[2], int mailbox);
  
  void		projectPosition(const Player *targ,const float t,float &x,float &y,float &z) const;
  void		getProjectedPosition(const Player *targ, float *projpos) const;
  
  private:
  const Player*	target;
  std::vector<RegionPoint>	path;
  int			pathIndex;
  float		timerForShot;
  bool		drivingForward;
  static std::vector<BzfRegion*>* obstacleList;
  
  /*variables added by Brendan McCabe*/
  Flag*               targetFlag;
  bool                flagsFound;
  bool                flocking;
  float               flock_direction[2];
  PlayerId            seeker;  //Id of the robot designated to capture flag
  MyPathfinder        pathfinder;
  std::vector<std::array<int, 3>>  my_path;
  bool                initialize_path;
  bool                need_path;
  float               pathfinder_mult;
  float               destination[3];
  bool                flag_change;
  const float* 	      endPoint;
 /*variables added by aidan akamine*/
  bool		      attacking;
  bool		      defending;
  int		      samePositionTimer;
  float		      getUnstuckPos[3];
  float		      moveToPosition[3];
};

#endif // BZF_ROBOT_PLAYER_H

// Local Variables: ***
// mode: C++ ***
// tab-width: 8 ***
// c-basic-offset: 2 ***
// indent-tabs-mode: t ***
// End: ***
// ex: shiftwidth=2 tabstop=8

