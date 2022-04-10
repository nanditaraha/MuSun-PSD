#include "FillHistBase.h"

FillHistBase::FillHistBase()
{
  PrintOutput = false;
}

FillHistBase::~FillHistBase()
{
}

/** Called on each module in succession after each tree entry
  * is read into the TMusunEvent in the passed MTAEvent.
  */
int FillHistBase::ProcessEntry(MTAEvent *mtaEvent)
{
  return 0;
}

/** Called before any entries are read from the input tree,
  * but after all FillHistBase modules have been created.
  */
int FillHistBase::Init(MTAEvent* mtaEvent)
{
  return 0;
}

/** Called after the last tree entry is read, but before
  * any modules are destructed or the input tree is closed.
  */
int FillHistBase::Terminate(MTAEvent* mtaEvent)
{
  return 0;
}
