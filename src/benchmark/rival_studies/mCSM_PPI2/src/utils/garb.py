#
#
# # protein = row["UniProt_ID"]
# # mutation = row["Mutation"]
# # interactor = row["Interactor_UniProt_ID"]
# # pdb_id = row["Template_cath_id_pdb"]
# # chain_id = row["Chain_id"]
#
# print(f"{protein=}")
# print(f"{mutation=}")
# print(f"{interactor=}")
# print(f"{pdb_id=}")
# print(f"{chain_id=}")
#
# upload_file(driver, pdb_id)
# # fill_pdb_accession(driver, pdb_id)
# fill_mutation(driver, mutation)
# fill_chain(driver, chain_id)
#
# submit_input(driver)
# process_validate_input_data(driver)
#
# url = get_url(driver)
#
# response = process_submission(driver)
# log.debug(f"Response: {response}")
# if response == ResponseMessages.COMPLETED:
#     log.info("Completed successfully.")
#     predicted_value, predicted_change = retrieve_completed_results(driver)
#     log.critical(f"predicted_value: {predicted_value}")
#     log.critical(f"predicted_change: {predicted_change}")
#     log.critical(f"url: {url}")
#
# else:
#     raise Exception
#
# time.sleep(100)
#
# # terminate_firefox_processes()
