#ifndef FillHistBase_h
#define FillHistBase_h 1

class MTAEvent;

/** The base class for MTA analysis modules.
  */
class FillHistBase
{
  public:
  FillHistBase();
  virtual ~FillHistBase();
  virtual int Init(MTAEvent *mtaEvent);
  virtual int ProcessEntry(MTAEvent *mtaEvent);
  virtual int Terminate(MTAEvent* mtaEvent);

  protected:
  bool PrintOutput;
};

#endif
