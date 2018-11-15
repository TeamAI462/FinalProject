/*
 * Defines the classes used for decision trees.
 *
 * Part of the Artificial Intelligence for Games system.
 *
 * Copyright (c) Ian Millington 2003-2006. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */
#include "dectree.h"

namespace aicore
{
  
  DecisionTreeNode* Decision::makeDecision(RobotPlayer* bot, float dt)
  {
    // Choose a branch based on the getBranch method
    if (getBranch(bot, dt)) {
      // Make sure its not null before recursing.
      if (trueBranch == NULL) return NULL;
      else return trueBranch->makeDecision(bot, dt);
    } else {
      // Make sure its not null before recursing.
      if (falseBranch == NULL) return NULL;
      else return falseBranch->makeDecision(bot, dt);
    }
  }
  
  DecisionTreeNode* DecisionPtr::makeDecision(RobotPlayer* bot, float dt)
  {
    // Choose a branch based on the getBranch method
    if ( getBranch(bot, dt) ) {
      // Make sure its not null before recursing.
      if (trueBranch == NULL) return NULL;
      else return trueBranch->makeDecision(bot, dt);
    } else {
      // Make sure its not null before recursing.
      if (falseBranch == NULL) return NULL;
      else return falseBranch->makeDecision(bot, dt);
    }
  }
  
  bool DecisionPtr::getBranch(RobotPlayer* bot, float dt)
  {
    return (bot->*decFuncPtr)(dt);
  }
  
  // Set up the trees
  void DecisionTrees::init()
  {
    /*start of lines modified by aidan akamine */
    /*******************doUpdateMotion(): evade or follow a*star search decision tree **********/
    /* amAlive: is robot alive?  */
    doUpdateMotionDecisions[0].decFuncPtr = &RobotPlayer::amAlive;
    doUpdateMotionDecisions[0].trueBranch = &doUpdateMotionDecisions[1];
    doUpdateMotionDecisions[0].falseBranch = &doUpdateMotionActions[0];
    
    /* isInDanger: is the robot about to get shot?*/
    doUpdateMotionDecisions[1].decFuncPtr = &RobotPlayer::isInDanger;
    doUpdateMotionDecisions[1].trueBranch = &doUpdateMotionActions[1];
    doUpdateMotionDecisions[1].falseBranch = &doUpdateMotionActions[2];
    
    
    /*******************doUpdate()  shooting decision tree  **********/
    /* amAlive: is robot alive?  */
    shootingDecisions[0].decFuncPtr = &RobotPlayer::amAlive;
    shootingDecisions[0].trueBranch = &shootingDecisions[1];
    shootingDecisions[0].falseBranch = &doUpdateMotionActions[0];
    
    /* readyToFire: is robot ready to fire a shot? */
    shootingDecisions[1].decFuncPtr = &RobotPlayer::readyToFire;
    shootingDecisions[1].trueBranch = &shootingDecisions[2];
    shootingDecisions[1].falseBranch = &doUpdateMotionActions[0];
    
    /* shotTimerReady: is the shot timer ready?*/
    shootingDecisions[2].decFuncPtr = &RobotPlayer::shotTimerReady;
    shootingDecisions[2].trueBranch = &shootingDecisions[3];
    shootingDecisions[2].falseBranch = &doUpdateMotionActions[0];
    
    /* shotWithinHalfTankLength if a shot is fired, will it be within half a tank length of enemy?*/
    shootingDecisions[3].decFuncPtr = &RobotPlayer::shotWithinHalfTankLenght;
    shootingDecisions[3].trueBranch = &shootingDecisions[4];
    shootingDecisions[3].falseBranch = &doUpdateMotionActions[0];
    
    /* isBuildingInBetween: is there a building blocking the shot?*/
    shootingDecisions[4].decFuncPtr = &RobotPlayer::isBuildingInBetween;
    shootingDecisions[4].trueBranch = &shootingDecisions[5];
    shootingDecisions[4].falseBranch = &doUpdateMotionActions[0];
    
    /* willShotHitTeammate: will a shot also hit a teammate? */
    shootingDecisions[5].decFuncPtr = &RobotPlayer::WillShotHitTeammate;
    shootingDecisions[5].trueBranch = &shootingActions[0];
    shootingDecisions[5].falseBranch = &shootingActions[1];
    
    /*********************doUpdate() decision tree on whether or not to drop flag *******/
    /* amAlive: is robot alive?  */
    holdingFlagDecisions[0].decFuncPtr = &RobotPlayer::amAlive;
    holdingFlagDecisions[0].trueBranch = &holdingFlagDecisions[1];
    holdingFlagDecisions[0].falseBranch = &doUpdateMotionActions[0];
    
    /*holdingFlag: is the tank holding a flag? */
    holdingFlagDecisions[1].decFuncPtr = &RobotPlayer::holdingFlag;
    holdingFlagDecisions[1].trueBranch = &holdingFlagDecisions[2];
    holdingFlagDecisions[1].falseBranch = &doUpdateMotionActions[0];
    
    /*isFlagSticky: is the flag the tank is holding sticky? */
    holdingFlagDecisions[2].decFuncPtr = &RobotPlayer::isFlagSticky;
    holdingFlagDecisions[2].trueBranch = &doUpdateMotionActions[0];
    holdingFlagDecisions[2].falseBranch = &holdingFlagDecisions[3];
    
    /*isFlagTeamFlag: is the flag the tank is holding a team flag? */
    holdingFlagDecisions[3].decFuncPtr = &RobotPlayer::isFlagTeamFlag;
    holdingFlagDecisions[3].trueBranch = &holdingFlagDecisions[4];
    holdingFlagDecisions[3].falseBranch = &holdingFlagActions[0];
    
    /*holdingMyTeamFlag: is the flag the tank is holding its own teamflag?*/
    holdingFlagDecisions[4].decFuncPtr = &RobotPlayer::holdingMyTeamFlag;
    holdingFlagDecisions[4].trueBranch = &holdingFlagActions[0];
    holdingFlagDecisions[4].falseBranch = &doUpdateMotionActions[0];
    
    
    /***** general actions  *********/
    doUpdateMotionActions[0].actFuncPtr = &RobotPlayer::doNothing;
    
    /*****actions for movement decision tree *********/
    doUpdateMotionActions[1].actFuncPtr = &RobotPlayer::evade;
    doUpdateMotionActions[2].actFuncPtr = &RobotPlayer::followPath;
    
    /*****actions for shooting decision tree *********/
    
    shootingActions[0].actFuncPtr = &RobotPlayer::Shooting;
    shootingActions[1].actFuncPtr = &RobotPlayer::teammateBlockingShot;
    
    
    /*****actions for flag dropping decision tree *********/
    holdingFlagActions[0].actFuncPtr = &RobotPlayer::dropFlag;
    
  }
  
  DecisionPtr DecisionTrees::doUpdateMotionDecisions[3];
  ActionPtr DecisionTrees::doUpdateMotionActions[4];
  
  DecisionPtr DecisionTrees::shootingDecisions[7];
  ActionPtr DecisionTrees::shootingActions[3];
  
  DecisionPtr DecisionTrees::holdingFlagDecisions[6];
  ActionPtr DecisionTrees::holdingFlagActions[2];
  
  /*end of lines modified by aidan akamine*/
}; // end of namespace
