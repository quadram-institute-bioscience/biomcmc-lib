#include <biomcmc.h> 
#include <check.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

#ifndef TEST_FILE_DIR
#define TEST_FILE_DIR "./files/"
#endif

char filename[2048] = TEST_FILE_DIR; // now we can memcpy() file names _after_ prefix_size
size_t prefix_size = strlen(TEST_FILE_DIR); // all modifications to filename[] come after prefix_size
newick_space nwk_spc;
speciestree sptre;
genetree gtre;

void
newick_space_setup_ortho_nwk (void)
{
  /* 4 small trees on same leafset; both unrooted and rooted, with branch lengths */
  memcpy(filename + prefix_size, "ortho.nwk", 10); 
  nwk_spc = new_newick_space_from_file (filename);
}

void
del_trees_teardown (void)
{
  del_newick_space (nwk_spc);
  del_speciestree (sptre);
  del_genetree (gtre);
}

START_TEST(new_single_topology_from_newick_file_function)
{
  topology t;
  memcpy(filename + prefix_size, "iqtree.nwk", 11); // 11 is sizeof file (or larger)
  t = new_single_topology_from_newick_file (filename);
  if (t->nleaves != 25) ck_abort_msg ("number of leaves disagree");
}
END_TEST

START_TEST(new_newick_space_from_file_ortho_nwk)
{
  if (nwk_spc->ntrees != 4) ck_abort_msg ("Problem reading 4 newick trees from ortho.nwk");
}
END_TEST

START_TEST(compare_ortho_nwk_unrooted_loop)
{
  int idx[5][3] = {{0,0,1},{0,1,1},{0,2,0},{2,3,1},{1,3,0}}; // [idx1, idx2, expected result]
  // [0] and [1] only differ by root location BUT topology needs same leaves AND same leaf order
  reorder_topology_leaves (nwk_spc->t[ idx[_i][0]]);
  reorder_topology_leaves (nwk_spc->t[ idx[_i][1]]);
  char_vector_link_address_if_identical (&(nwk_spc->t[ idx[_i][0]]->taxlabel), &(nwk_spc->t[ idx[_i][1]]->taxlabel));
  int res = topology_is_equal_unrooted (nwk_spc->t[ idx[_i][0]], nwk_spc->t[ idx[_i][1]], false);
  if (res != idx[_i][2]) ck_abort_msg ("unrooted comparison with ortho.nwk");
}
END_TEST

START_TEST(compare_ortho_nwk_rooted_loop)
{
  int idx[5][2] = {{0,1},{0,2},{0,3},{1,2},{1,3}}; // all pairs should be distinct except {2,3} by chance of remove_multifurcation
  // [0] and [1] only differ by root location BUT topology needs same leaves AND same leaf order
  reorder_topology_leaves (nwk_spc->t[ idx[_i][0]]);
  reorder_topology_leaves (nwk_spc->t[ idx[_i][1]]);
  char_vector_link_address_if_identical (&(nwk_spc->t[ idx[_i][0]]->taxlabel), &(nwk_spc->t[ idx[_i][1]]->taxlabel));
  int res = topology_is_equal (nwk_spc->t[ idx[_i][0]], nwk_spc->t[ idx[_i][1]]);
  if (res) ck_abort_msg ("rooted comparison with ortho.nwk");
}
END_TEST

START_TEST(new_speciestree_function)
{
  int idx[3][6] = {{1,0,0,0,0,0}, {2,13,106,80,44,103}, {3,13,106,80,44,103}}; // no SPR since not exact
  bool is_correct = true;
  sptre = new_speciestree (nwk_spc->t[0], NULL);
  gtre = new_genetree (nwk_spc->t[idx[_i][0]], sptre);
  genetree_reconcile_speciestree (gtre, sptre); 
  genetree_dSPR_speciestree (gtre, sptre, 2);
  is_correct = is_correct && (gtre->rec->ndups == idx[_i][1]);
  is_correct = is_correct && (gtre->rec->nloss == idx[_i][2]);
  is_correct = is_correct && (gtre->rec->ndcos == idx[_i][3]);
  is_correct = is_correct && (gtre->split->rf    == idx[_i][4]);
  is_correct = is_correct && (gtre->split->hdist == idx[_i][5]);
  if (!is_correct) ck_abort_msg ("Reconciliation/splitset distances distinct from expected %d", gtre->rec->ndups);
}
END_TEST

Suite * topology_suite(void)
{
  Suite *s;
  TCase *tc_case;

  s = suite_create("Topology");

  tc_case = tcase_create("read_trees");
  tcase_add_test(tc_case, new_single_topology_from_newick_file_function);
  suite_add_tcase(s, tc_case);
  tc_case = tcase_create("read_newick_space");
  tcase_add_checked_fixture(tc_case, newick_space_setup_ortho_nwk, del_trees_teardown); // unchecked -> once per case; checked -> per unit
  tcase_add_test(tc_case, new_newick_space_from_file_ortho_nwk);
  tcase_add_loop_test (tc_case, compare_ortho_nwk_unrooted_loop, 0, 5); // loops, using index _i
  tcase_add_loop_test (tc_case, compare_ortho_nwk_rooted_loop, 0, 5); 
  suite_add_tcase(s, tc_case);
  tc_case = tcase_create("create_gene_species_trees");
  tcase_add_checked_fixture(tc_case, newick_space_setup_ortho_nwk, del_trees_teardown); // unchecked -> once per case; checked -> per unit
  tcase_add_loop_test (tc_case, new_speciestree_function, 0, 3); 
  suite_add_tcase(s, tc_case);

  return s;
}

int main(void)
{
  int number_failed;
  SRunner *sr;

  sr = srunner_create (topology_suite());
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed > 0) ? TEST_FAILURE:TEST_SUCCESS;
}
