
/**
 * A very simple unit testing class.
 *
 * @todo have each test() forked into a different
 *  process so that seg faults can be caught and the
 *  tests continued past the offending (and consequently
 *  marked) code.
 * 
 * Copyright (C) 2002, 2003 David M. Doolin
 *
 * This file is part of Geotechnica.
 *
 * Geotechnica is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * Geotechnica is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
 * for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with Geotechnica; see the file COPYING.  If not, write to the Free
 * Software Foundation, 59 Temple Place - Suite 330, Boston,
 * MA  02111-1307, USA.
 */


#ifndef __GEO_UNITTEST_H__
#define __GEO_UNITTEST_H__


/** Each unit test builds a table that can be 
 * traversed to exercise each function in the 
 * component.
 */
typedef struct _testfunc {

  bool       (*test)(void);
  const char * testname;
} TestFunc;



class UnitTest {

 public:

  void register_test_functions (TestFunc * testfunc);

 /**  Perform all of the unit testing specified in a 
  *  table.
  *
  *  @param TestFunc * an array of test functions that 
  *  are invoked sequentially from the calling function.
  *
  *  @return int TRUE if everything passes the unit test,
  *  false if any function fails its unit test.
  */
  bool test                    ();

  void print_header            (void * stream, 
                                const char * tag);

 private:

  TestFunc * testfunc;

};



/** 
 * Individual test functions can be prototyped 
 * here for inclusion into other test program drivers.
 */


#endif  /* __GEO_UNITTEST_H__ */
