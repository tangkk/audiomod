#ifndef _DATATYPE_H
#define _DATATYPE_H

struct NoteStru
{
	int			iStart;
	int			iDuration;
	int			iNote;
} ;

struct fNoteStru
{
	int			iStart;
	int			iDuration;
	float		iNote;	// in Hz
} ;

struct fNoteStru0
{
	float			iStart;
	float			iDuration;
	float		iNote;	// in Hz
} ;

struct MidGenStru
{
	int			iTrack;
	int			iInstrument;
	int			iPan;
	float		iVolume; //range = [0-1]
	int			iTempo;
} ;

/*
struct NoteGenApp1Stru
{
	vector< pair<float, float> > smoothPitch;
	vector< float > mLevel;
} ;
*/

#endif /* _DATATYPE_H */